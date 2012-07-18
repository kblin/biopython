# Copyright 2012 by Kai Blin.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.SearchIO parser for HMMER 2 text output."""

import re
from Bio.SearchIO._objects import QueryResult, Hit, HSP

_HSP_ALIGN_LINE = re.compile(r'(\S+):\s+domain (\d+) of (\d+)')

class Hmmer2TextIterator(object):
    """Iterator for the HMMER 2.0 text output."""

    def __init__(self, handle):
        self.handle = handle
        self.buf = []
        self._meta = self.parse_preamble()

    def __iter__(self):
        for qresult in self.parse_qresult():
            qresult.program = self._meta.get('program')
            qresult.target = self._meta.get('target')
            qresult.version = self._meta.get('version')
            yield qresult

    def read_next(self):
        """Return the next non-empty line, trailing whitespace removed"""
        if len(self.buf) > 0:
            return self.buf.pop()
        self.line = self.handle.readline()
        while self.line and not self.line.strip():
            self.line = self.handle.readline()
        if self.line:
            self.line = self.line.rstrip()
        return self.line

    def push_back(self, line):
        """Un-read a line that should not be parsed yet"""
        self.buf.append(line)

    def parse_key_value(self):
        """Parse key-value pair separated by colon (:)"""
        key, value = self.line.split(':')
        return key.strip(), value.strip()

    def parse_preamble(self):
        """Parse HMMER2 preamble."""
        meta = {}
        state = "GENERIC"
        while self.read_next():
            if state == "GENERIC":
                if self.line.startswith('hmm'):
                    meta['program'] = self.line.split('-')[0].strip()
                elif self.line.startswith('HMMER is'):
                    continue
                elif self.line.startswith('HMMER'):
                    meta['version'] = self.line.split()[1]
                elif self.line.count('-') == 36:
                    state = "OPTIONS"
                continue

            assert state == "OPTIONS"
            assert 'program' in meta

            if self.line.count('-') == 32:
                break

            key, value = self.parse_key_value()
            if meta['program'] == 'hmmsearch':
                if key == 'Sequence database':
                    meta['target'] = value
                    continue
            elif meta['program'] == 'hmmpfam':
                if key == 'HMM file':
                    meta['target'] = value
                    continue
            meta[key] = value

        return meta

    def parse_qresult(self):
        """Parse a HMMER2 query block."""
        while self.read_next():
            if not self.line.startswith('Query'):
                raise StopIteration()
            _, id_ = self.parse_key_value()
            self.qresult = QueryResult(id_)

            while self.read_next() and not self.line.startswith('Scores'):
                if self.line.startswith('Accession'):
                    self.qresult.acc = self.parse_key_value()[1]
                if self.line.startswith('Description'):
                    self.qresult.desc = self.parse_key_value()[1]

            self.parse_hits()
            self.parse_hsps()
            self.parse_hsp_alignments()

            while self.read_next() and self.line != '//':
                pass

            yield self.qresult


    def parse_hits(self):
        """Parse a HMMER2 hit block, beginning with the hit table."""

        while self.read_next():
            if self.line.startswith('Parsed'):
                break

            if self.line.startswith('Sequence') or \
               self.line.startswith('Model') or \
               self.line.startswith('-------- '):
                continue

            fields = self.line.split()
            id_ = fields.pop(0)
            n = int(fields.pop())
            evalue = float(fields.pop())
            score = float(fields.pop())
            desc = ' '.join(fields)


            hit = Hit(id_, self.qresult.id)
            hit.evalue = evalue
            hit.score = score
            hit.desc = desc
            hit.n = n
            self.qresult.append(hit)


    def parse_hsps(self):
        """Parse a HMMER2 hsp block, beginning with the hsp table."""
        while self.read_next():
            if self.line.startswith('Alignments') or \
               self.line.startswith('Histogram') or \
               self.line == '//':
                break
            if self.line.startswith('Model') or \
               self.line.startswith('Sequence') or \
               self.line.startswith('--------'):
                continue

            id_, domain, seq_f, seq_t, seq_compl, hmm_f, hmm_t, hmm_compl, \
            score, evalue = self.line.split()

            hsp = HSP(id_, self.qresult.id)
            hsp.evalue = float(evalue)
            hsp.bitscore = float(score)
            hsp.domain_index = int(domain.split('/')[0])
            if self._meta['program'] == 'hmmpfam':
                hsp.hit_start = int(hmm_f)
                hsp.hit_end = int(hmm_t)
                hsp.query_start = int(seq_f)
                hsp.query_end = int(seq_t)
            elif self._meta['program'] == 'hmmsearch':
                hsp.query_start = int(seq_f)
                hsp.query_end = int(seq_t)
                hsp.hit_start = int(hmm_f)
                hsp.hit_end = int(hmm_t)


            self.qresult[id_].append(hsp)



    def parse_hsp_alignments(self):
        """Parse a HMMER2 HSP alignment block."""
        if not self.line.startswith('Alignments'):
            return

        while self.read_next():
            if self.line == '//' or self.line.startswith('Histogram'):
                break

            match = re.search(_HSP_ALIGN_LINE, self.line)
            if match is None:
                continue

            id_ = match.group(1)
            idx = int(match.group(2))
            n = int(match.group(3))

            hit = self.qresult[id_]
            if hit.n != n:
                continue

            hsp = hit[idx-1]

            hmmseq = ''
            consensus = ''
            otherseq = ''
            while self.read_next() and self.line.startswith(' '):
                if self.line[19] == '*':
                    hmmseq += self.line[22:]
                else:
                    hmmseq += self.line[19:]

                if not self.read_next():
                    break
                consensus += self.line[19:].strip()

                if not self.read_next():
                    break
                otherseq += self.line[19:].split()[0].strip()

            self.push_back(self.line)

            if self._meta['program'] == 'hmmpfam':
                hsp.hit = hmmseq
                hsp.query = otherseq
            else:
                hsp.hit = otherseq
                hsp.query = hmmseq
