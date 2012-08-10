# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER domain table output format."""

from itertools import chain

from Bio.SearchIO._objects import QueryResult, Hit, HSP, HSPFragment
from hmmertab import HmmerTabIterator, HmmerTabIndexer


class HmmerDomtabIterator(HmmerTabIterator):

    """Base hmmer-domtab iterator."""

    def parse_result_row(self):
        """Returns a dictionary of parsed row values."""
        assert self.line
        cols = filter(None, self.line.strip().split(' '))
        # if len(cols) > 23, we have extra description columns
        # combine them all into one string in the 19th column
        if len(cols) > 23:
            cols[22] = ' '.join(cols[22:])
        elif len(cols) < 23:
            cols.append('')
            assert len(cols) == 23

        # assign parsed column data into qresult, hit, and hsp dicts
        qresult = {}
        qresult['id'] = cols[3]                 # query name
        qresult['acc'] = cols[4]                # query accession
        qresult['seq_len'] = int(cols[5])       # qlen
        hit = {}
        hit['id'] = cols[0]                     # target name
        hit['acc'] = cols[1]                    # target accession
        hit['seq_len'] = int(cols[2])           # tlen
        hit['evalue'] = float(cols[6])          # evalue
        hit['bitscore'] = float(cols[7])        # score
        hit['bias'] = float(cols[8])            # bias
        hit['description'] = cols[22]           # description of target
        hsp = {}
        hsp['domain_index'] = int(cols[9])      # # (domain number)
        # not parsing cols[10] since it's basically len(hit)
        hsp['evalue_cond'] = float(cols[11])    # c-evalue
        hsp['evalue'] = float(cols[12])         # i-evalue
        hsp['bitscore'] = float(cols[13])       # score
        hsp['bias'] = float(cols[14])           # bias
        hsp['env_start'] = int(cols[19]) - 1    # env from
        hsp['env_end'] = int(cols[20])          # env to
        hsp['acc_avg'] = float(cols[21])        # acc
        frag = {}
        # strand is always 0, since HMMER now only handles protein
        frag['hit_strand'] = frag['query_strand'] = 0
        frag['hit_start'] = int(cols[15]) - 1    # hmm from
        frag['hit_end'] = int(cols[16])          # hmm to
        frag['query_start'] = int(cols[17]) - 1  # ali from
        frag['query_end'] = int(cols[18])        # ali to

        # switch hmm<-->ali coordinates if hmm is not hit
        if not self.hmm_as_hit:
            frag['hit_end'], frag['query_end'] = \
                    frag['query_end'], frag['hit_end']
            frag['hit_start'], frag['query_start'] = \
                    frag['query_start'], frag['hit_start']

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp, 'frag': frag}

    def parse_qresult(self):
        """Generator function that returns QueryResult objects."""
        # state values, determines what to do for each line
        state_EOF = 0
        state_QRES_NEW = 1
        state_QRES_SAME = 3
        state_HIT_NEW = 2
        state_HIT_SAME = 4
        # dummies for initial states
        qres_state = None
        hit_state = None
        file_state = None
        # dummies for initial id caches
        prev_qid = None
        prev_hid = None
        # dummies for initial parsed value containers
        cur, prev = None, None
        hit_list, hsp_list = [], []

        while True:
            # store previous line's parsed values, for every line after the 1st
            if cur is not None:
                prev = cur
                prev_qid = cur_qid
                prev_hid = cur_hid
            # only parse the line if it's not EOF
            if self.line:
                cur = self.parse_result_row()
                cur_qid = cur['qresult']['id']
                cur_hid = cur['hit']['id']
            else:
                file_state = state_EOF
                # mock ID values since the line is empty
                cur_qid, cur_hid = None, None

            # get the state of hit and qresult
            if prev_qid != cur_qid:
                qres_state = state_QRES_NEW
            else:
                qres_state = state_QRES_SAME
            # new hits are hits with different ids or hits in a new qresult
            if prev_hid != cur_hid or qres_state == state_QRES_NEW:
                hit_state = state_HIT_NEW
            else:
                hit_state = state_HIT_SAME

            # start creating objects after the first line (i.e. prev is filled)
            if prev is not None:
                # each line is basically an HSP with one HSPFragment
                frag = HSPFragment(prev_hid, prev_qid)
                for attr, value in prev['frag'].items():
                    setattr(frag, attr, value)
                hsp = HSP([frag])
                for attr, value in prev['hsp'].items():
                    setattr(hsp, attr, value)
                hsp_list.append(hsp)

                # create hit object when we've finished parsing all its hsps
                # i.e. when hit state is state_HIT_NEW
                if hit_state == state_HIT_NEW:
                    hit = Hit(prev_hid, prev_qid, hsps=hsp_list)
                    for attr, value in prev['hit'].items():
                        setattr(hit, attr, value)
                    hit_list.append(hit)
                    hsp_list = []

                # create qresult and yield if we're at a new qresult or EOF
                if qres_state == state_QRES_NEW or file_state == state_EOF:
                    qresult = QueryResult(prev_qid, hits=hit_list)
                    for attr, value in prev['qresult'].items():
                        setattr(qresult, attr, value)
                    yield qresult
                    # if current line is EOF, break
                    if file_state == state_EOF:
                        break
                    hit_list = []

            self.line = self.handle.readline()


class HmmerDomtabHmmhitIterator(HmmerDomtabIterator):

    """Parser for the HMMER domain table format that assumes HMM profile
    coordinates are hit coordinates."""

    hmm_as_hit = True


class HmmerDomtabHmmqueryIterator(HmmerDomtabIterator):

    """Parser for the HMMER domain table format that assumes HMM profile
    coordinates are query coordinates."""

    hmm_as_hit = False


class HmmerDomtabHmmhitIndexer(HmmerTabIndexer):

    """Indexer class for HMMER domain table output that assumes HMM profile
    coordinates are hit coordinates."""

    _parser = HmmerDomtabHmmhitIterator
    _query_id_idx = 3


class HmmerDomtabHmmqueryIndexer(HmmerTabIndexer):

    """Indexer class for HMMER domain table output that assumes HMM profile
    coordinates are query coordinates."""

    _parser = HmmerDomtabHmmqueryIterator
    _query_id_idx = 3


class HmmerDomtabHmmhitWriter(object):

    """Writer for hmmer-domtab output format which writes hit coordinates
    as HMM profile coordinates."""

    hmm_as_hit = True

    def __init__(self, handle):
        self.handle = handle

    def write_file(self, qresults):
        """Writes to the handle.

        Returns a tuple of how many QueryResult, Hit, and HSP objects were written.

        """
        handle = self.handle
        qresult_counter, hit_counter, hsp_counter = 0, 0, 0

        try:
            first_qresult = qresults.next()
        except StopIteration:
            handle.write(self.build_header())
        else:
            # write header
            handle.write(self.build_header(first_qresult))
            # and then the qresults
            for qresult in chain([first_qresult], qresults):
                if qresult:
                    handle.write(self.build_row(qresult))
                    qresult_counter += 1
                    hit_counter += len(qresult)
                    hsp_counter += sum([len(hit) for hit in qresult])

        return qresult_counter, hit_counter, hsp_counter

    def build_header(self, first_qresult=None):
        """Returns the header string of a domain HMMER table output."""

        # calculate whitespace required
        # adapted from HMMER's source: src/p7_tophits.c#L1157
        if first_qresult:
            #qnamew = max(20, len(first_qresult.id))
            qnamew = 20
            tnamew = max(20, len(first_qresult[0].id))
            try:
                qaccw = max(10, len(first_qresult.acc))
                taccw = max(10, len(first_qresult[0].acc))
            except AttributeError:
                qaccw, taccw = 10, 10
        else:
            qnamew, tnamew, qaccw, taccw = 20, 20, 10, 10

        header = ''
        header += "#%*s %22s %40s %11s %11s %11s\n" % \
                (tnamew+qnamew-1+15+taccw+qaccw, "", "--- full sequence ---", \
                "-------------- this domain -------------", "hmm coord", \
                "ali coord", "env coord")
        header += "#%-*s %-*s %5s %-*s %-*s %5s %9s %6s %5s %3s %3s %9s " \
                "%9s %6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n" % (tnamew-1, \
                " target name", taccw, "accession", "tlen", qnamew, \
                "query name", qaccw, "accession", "qlen", "E-value", "score", \
                "bias", "#", "of", "c-Evalue", "i-Evalue", "score", "bias", \
                "from", "to", "from", "to", "from", "to", "acc", \
                "description of target")
        header += "#%*s %*s %5s %*s %*s %5s %9s %6s %5s %3s %3s %9s %9s " \
                "%6s %5s %5s %5s %5s %5s %5s %5s %4s %s\n" % (tnamew-1, \
                "-------------------", taccw, "----------", "-----", \
                qnamew, "--------------------", qaccw, "----------", \
                "-----", "---------", "------", "-----", "---", "---", \
                "---------", "---------", "------", "-----", "-----", "-----", \
                "-----", "-----", "-----", "-----", "----", \
                "---------------------")

        return header

    def build_row(self, qresult):
        """Returns a string or one row or more of the QueryResult object."""
        rows = ''

        # calculate whitespace required
        # adapted from HMMER's source: src/p7_tophits.c#L1083
        qnamew = max(20, len(qresult.id))
        tnamew = max(20, len(qresult[0].id))
        try:
            qaccw = max(10, len(qresult.acc))
            taccw = max(10, len(qresult[0].acc))
            qresult_acc = qresult.acc
        except AttributeError:
            qaccw, taccw = 10, 10
            qresult_acc = '-'

        for hit in qresult:
            
            # try to get hit accession
            try:
                hit_acc = hit.acc
            except AttributeError:
                hit_acc = '-'

            for hsp in hit.hsps:
                if self.hmm_as_hit:
                    hmm_to = hsp.hit_end
                    hmm_from = hsp.hit_start + 1
                    ali_to = hsp.query_end
                    ali_from = hsp.query_start + 1
                else:
                    hmm_to = hsp.query_end
                    hmm_from = hsp.query_start + 1
                    ali_to = hsp.hit_end
                    ali_from = hsp.hit_start + 1

                rows += "%-*s %-*s %5d %-*s %-*s %5d %9.2g %6.1f %5.1f %3d %3d" \
                " %9.2g %9.2g %6.1f %5.1f %5d %5d %5ld %5ld %5d %5d %4.2f %s\n" % \
                (tnamew, hit.id, taccw, hit_acc, hit.seq_len, qnamew, qresult.id, \
                qaccw, qresult_acc, qresult.seq_len, hit.evalue, hit.bitscore, \
                hit.bias, hsp.domain_index, len(hit.hsps), hsp.evalue_cond, hsp.evalue, \
                hsp.bitscore, hsp.bias, hmm_from, hmm_to, ali_from, ali_to, \
                hsp.env_start + 1, hsp.env_end, hsp.acc_avg, hit.description)

        return rows


class HmmerDomtabHmmqueryWriter(HmmerDomtabHmmhitWriter):

    """Writer for hmmer-domtab output format which writes query coordinates
    as HMM profile coordinates."""

    hmm_as_hit = False


def _test():
    """Run the Bio.SearchIO.HmmerIO.hmmerdomtab module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os

    test_dir = 'Tests'

    if os.path.isdir(os.path.join('..', '..', test_dir)):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join('..', '..', test_dir))
        doctest.testmod()
        os.chdir(cur_dir)
        print "Done"


# if not used as a module, run the doctest
if __name__ == "__main__":
    _test()
