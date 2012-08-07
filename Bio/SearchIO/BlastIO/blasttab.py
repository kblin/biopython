# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ tab output format, with or without comments."""

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._objects import QueryResult, Hit, HSP, HSPFragment


__all__ = ['BlastTabIndexer', 'BlastTabIterator', 'BlastTabWriter']

# longname-shortname map
# maps the column names shown in a commented output to its short name
# (the one used in the command line)
_LONG_SHORT_MAP = {
    'query id': 'qseqid',
    'query acc.': 'qacc',
    'query acc.ver': 'qaccver',
    'query length': 'qlen',
    'subject id': 'sseqid',
    'subject acc.': 'sacc',
    'subject acc.ver': 'saccver',
    'subject length': 'slen',
    'alignment length': 'length',
    'bit score': 'bitscore',
    'score': 'score',
    'evalue': 'evalue',
    'identical': 'nident',
    '% identity': 'pident',
    'positives': 'positive',
    '% positives': 'ppos',
    'mismatches': 'mismatch',
    'gaps': 'gaps',
    'q. start': 'qstart',
    'q. end': 'qend',
    's. start': 'sstart',
    's. end': 'send',
    'query frame': 'qframe',
    'sbjct frame': 'sframe',
    'query/sbjct frames': 'frames',
    'query seq': 'qseq',
    'subject seq': 'sseq',
    'gap opens': 'gapopen',
    'query gi': 'qgi',
    'subject ids': 'sallseqid',
    'subject gi': 'sgi',
    'subject gis': 'sallgi',
    'BTOP': 'btop',
}

# column to class attribute map
_COLUMN_QRESULT = {
    'qseqid': ('id', str),
    'qacc': ('acc', str),
    'qaccver': ('acc_ver', str),
    'qlen': ('seq_len', int),
    'qgi': ('gi', str),
}
_COLUMN_HIT = {
    'sseqid': ('id', str),
    'sallseqid': ('id_all', str),
    'sacc': ('acc', str),
    'saccver': ('acc_ver', str),
    'sgi': ('gi', str),
    'sallgi': ('gi_all', str),
    'slen': ('seq_len', int),
}
_COLUMN_HSP = {
    'bitscore': ('bitscore', float),
    'score': ('bitscore_raw', int),
    'evalue': ('evalue', float),
    'nident': ('ident_num', int),
    'pident': ('ident_pct', float),
    'positive': ('pos_num', int),
    'ppos': ('pos_pct', float),
    'mismatch': ('mismatch_num', int),
    'gaps': ('gap_num', int),
    'gapopen': ('gapopen_num', int),
    'btop': ('btop', str),
}
_COLUMN_FRAG = {
    'length': ('aln_span', int),
    'qstart': ('query_start', int),
    'qend': ('query_end', int),
    'sstart': ('hit_start', int),
    'send': ('hit_end', int),
    'qframe': ('query_frame', int),
    'sframe': ('hit_frame', int),
    'frames': ('frames', str),
    'qseq': ('query', str),
    'sseq': ('hit', str),
}
_SUPPORTED_FIELDS = set(_COLUMN_QRESULT.keys() + _COLUMN_HIT.keys() + \
        _COLUMN_HSP.keys() + _COLUMN_FRAG.keys())

# column order in the non-commented tabular output variant
# values must be keys inside the column-attribute maps above
_DEFAULT_FIELDS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', \
        'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
# one field from each of the following sets must exist in order for the
# parser to work
_MIN_QUERY_FIELDS = set(['qseqid', 'qacc', 'qaccver'])
_MIN_HIT_FIELDS = set(['sseqid', 'sacc', 'saccver'])


class BlastTabIterator(object):

    """Main parser for the BLAST tabular format."""

    def __init__(self, handle, comments=False, fields=_DEFAULT_FIELDS):
        self.handle = handle
        self.has_comments = comments
        self.fields = self._prep_fields(fields)
        self.line = self.handle.readline().strip()

    def __iter__(self):
        # stop iteration if file has no lines
        if not self.line:
            raise StopIteration
        # determine which iterator to use
        elif self.has_comments:
            iterfunc = self._parse_qresult_with_comments
        else:
            iterfunc = self._parse_qresult

        for qresult in iterfunc():
            yield qresult

    def _prep_fields(self, fields):
        """Validates and formats the given fields for use by the parser."""
        # cast into list if fields is a space-separated string
        if isinstance(fields, basestring):
            fields = fields.strip().split(' ')
        # blast allows 'std' as a proxy for the standard default lists
        # we want to transform 'std' to its proper column names
        if 'std' in fields:
            idx = fields.index('std')
            fields = fields[:idx] + _DEFAULT_FIELDS + fields[idx+1:]
        # if set(fields) has a null intersection with minimum required
        # fields for hit and query, raise an exception
        if not set(fields).intersection(_MIN_QUERY_FIELDS) or \
                not set(fields).intersection(_MIN_HIT_FIELDS):
            raise ValueError("Required query and/or hit ID field not found.")

        return fields

    def _parse_qresult_with_comments(self):
        """Iterator returning `QueryResult` objects from a commented file."""
        while True:
            comments = self._parse_comments()
            if comments:
                try:
                    self.fields = comments['fields']
                    empty_qres = False
                except KeyError:
                    # no fields means the query has no results
                    assert 'fields' not in comments
                    empty_qres = True

                if empty_qres:
                    qresult = QueryResult('')
                    for key, value in comments.items():
                        setattr(qresult, key, value)
                    yield qresult
                else:
                    for qresult in self._parse_qresult():
                        for key, value in comments.items():
                            setattr(qresult, key, value)
                        yield qresult
            else: break

    def _parse_comments(self):
        """Returns a dictionary containing tab file comments."""
        comments = {}
        while True:
            # parse program and version
            # example: # BLASTX 2.2.26+
            if 'BLAST' in self.line and 'processed' not in self.line:
                program_line = self.line[len(' #'):].split(' ')
                comments['program'] = program_line[0].lower()
                comments['version'] = program_line[1]
            # parse query id and description (if available)
            # example: # Query: gi|356995852 Mus musculus POU domain
            elif 'Query' in self.line:
                query_line = self.line[len('# Query: '):].split(' ', 1)
                comments['id'] = query_line[0]
                if len(query_line) == 2:
                    comments['description'] = query_line[1]
            # parse target database
            # example: # Database: db/minirefseq_protein
            elif 'Database' in self.line:
                comments['target'] = self.line[len('# Database: '):]
            # parse RID (from remote searches)
            elif 'RID' in self.line:
                comments['rid'] = self.line[len('# RID: '):]
            # parse column order, required for parsing the result lines
            # example: # Fields: query id, query gi, query acc., query length
            elif 'Fields' in self.line:
                comments['fields'] = self._parse_fields_line()
            # if the line has these strings, it's either the end of a comment
            # or the end of a file, so we return all the comments we've parsed
            elif ' hits found' in self.line or 'processed' in self.line:
                self.line = self.handle.readline().strip()
                return comments

            self.line = self.handle.readline()

            if not self.line:
                return
            else:
                self.line = self.line.strip()

    def _parse_fields_line(self):
        """Returns a list of column short names from the 'Fields'
        comment line."""
        raw_field_str = self.line[len('# Fields: '):]
        long_fields = raw_field_str.split(', ')
        fields = [_LONG_SHORT_MAP[long_name] for long_name in long_fields]
        return self._prep_fields(fields)

    def _parse_result_row(self):
        """Returns a dictionary of parsed row values."""
        fields = self.fields
        columns = self.line.strip().split('\t')
        assert len(fields) == len(columns), "Expected %i columns, found: " \
            "%i" % (len(fields), len(columns))

        qresult, hit, hsp, frag = {}, {}, {}, {}
        for idx, value in enumerate(columns):
            sname = fields[idx]
            # flag to check if any of the _COLUMNs contain sname
            in_mapping = False
            # iterate over each dict, mapping pair to determine
            # attribute name and value of each column
            for parsed_dict, mapping in (
                    (qresult, _COLUMN_QRESULT),
                    (hit, _COLUMN_HIT),
                    (hsp, _COLUMN_HSP),
                    (frag, _COLUMN_FRAG)):
                # process parsed value according to mapping
                if sname in mapping:
                    attr_name, caster = mapping[sname]
                    if caster is not str:
                        value = caster(value)
                    parsed_dict[attr_name] = value
                    in_mapping = True
            # make sure that any unhandled field is not supported
            if not in_mapping:
                assert sname not in _SUPPORTED_FIELDS

        return {'qresult': qresult, 'hit': hit, 'hsp': hsp, 'frag': frag}

    def _get_id(self, parsed):
        """Returns the value used for a QueryResult or Hit ID from a parsed row."""
        # use 'id', with 'acc' and 'acc_ver' fallbacks
        # one of these must have a value since we've checked whether
        # they exist or not when parsing the comments
        id_cache = parsed.get('id')
        if id_cache is None:
            id_cache = parsed.get('acc')
        if id_cache is None:
            id_cache = parsed.get('acc_ver')

        return id_cache

    def _parse_qresult(self):
        """Generator function that returns QueryResult objects."""
        # state values, used to determine what to do with each line
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
            # store previous line's parsed values if we've past the first line
            if cur is not None:
                prev = cur
                prev_qid = cur_qid
                prev_hid = cur_hid
            # only parse the line if it's not EOF or not a comment line
            if self.line and not self.line.startswith('#'):
                cur = self._parse_result_row()
                cur_qid = self._get_id(cur['qresult'])
                cur_hid = self._get_id(cur['hit'])
            else:
                file_state = state_EOF
                # mock values for cur_qid and cur_hid since the line is empty
                cur_qid, cur_hid = None, None

            # get the state of hit and qresult
            if prev_qid != cur_qid:
                qres_state = state_QRES_NEW
            else:
                qres_state = state_QRES_SAME
            # new hits are hits with different id or hits in a new qresult
            if prev_hid != cur_hid or qres_state == state_QRES_NEW:
                hit_state = state_HIT_NEW
            else:
                hit_state = state_HIT_SAME

            # we're creating objects for the previously parsed line(s),
            # so nothing is done in the first parsed line (prev == None)
            if prev is not None:
                # every line is essentially an HSP with one fragment, so we
                # create both of these for every line
                frag = HSPFragment(prev_hid, prev_qid)
                for attr, value in prev['frag'].items():
                    # adjust coordinates to Python range
                    # NOTE: this requires both start and end coords to be
                    # present, otherwise a KeyError will be raised.
                    # Without this limitation, we might misleadingly set the
                    # start / end coords
                    for seq_type in ('query', 'hit'):
                        if attr == seq_type + '_start':
                            value = min(value,
                                    prev['frag'][seq_type + '_end']) - 1
                        elif attr == seq_type + '_end':
                            value = max(value,
                                    prev['frag'][seq_type + '_start'])
                    setattr(frag, attr, value)
                # strand and frame setattr require the full parsed values
                # to be set first
                for seq_type in ('hit', 'query'):
                    # try to set hit and query frame
                    frame = self._get_frag_frame(frag, seq_type,
                            prev['frag'])
                    setattr(frag, '%s_frame' % seq_type, frame)
                    # try to set hit and query strand
                    strand = self._get_frag_strand(frag, seq_type,
                            prev['frag'])
                    setattr(frag, '%s_strand' % seq_type, strand)

                hsp = HSP([frag])
                for attr, value in prev['hsp'].items():
                    setattr(hsp, attr, value)
                hsp_list.append(hsp)

                # create hit and append to temp hit container if hit_state
                # says we're not at the same hit or at a new query
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

            self.line = self.handle.readline().strip()

    def _get_frag_frame(self, frag, seq_type, parsedict):
        """Returns `HSPFragment` frame given the object, its sequence type,
        and its parsed dictionary values."""
        assert seq_type in ('query', 'hit')
        frame = getattr(frag, '%s_frame' % seq_type, None)
        if frame is not None:
            return frame
        else:
            if 'frames' in parsedict:
                # frames is 'x1/x2' string, x1 is query frame, x2 is hit frame
                idx = 0 if seq_type == 'query' else 1
                return int(parsedict['frames'].split('/')[idx])
            # else implicit None return

    def _get_frag_strand(self, frag, seq_type, parsedict):
        """Returns `HSPFragment` strand given the object, its sequence type,
        and its parsed dictionary values."""
        # NOTE: this will never set the strands as 0 for protein
        # queries / hits, since we can't detect the blast flavors
        # from the columns alone.
        assert seq_type in ('query', 'hit')
        strand = getattr(frag, '%s_strand' % seq_type, None)
        if strand is not None:
            return strand
        else:
            # using parsedict instead of the fragment object since
            # we need the unadjusted coordinated values
            start = parsedict.get('%s_start' % seq_type)
            end = parsedict.get('%s_end' % seq_type)
            if start is not None and end is not None:
                return 1 if start <= end else -1
            # else implicit None return


class BlastTabIndexer(SearchIndexer):

    """Indexer class for BLAST+ tab output."""

    _parser = BlastTabIterator

    def __init__(self, filename, comments=False, fields=_DEFAULT_FIELDS):
        SearchIndexer.__init__(self, filename, comments=comments, fields=fields)
        self._has_comments = comments
        self._handle.seek(0)

        # if the file doesn't have comments,
        # get index of column used as the key (qseqid / qacc / qaccver)
        if not self._has_comments:
            if 'qseqid' in fields:
                self._key_idx = fields.index('qseqid')
            elif 'qacc' in fields:
                self._key_idx = fields.index('qacc')
            elif 'qaccver' in fields:
                self._key_idx = fields.index('qaccver')
            else:
                raise ValueError("Custom fields is missing an ID column. "
                        "One of these must be present: 'qseqid', 'qacc', or 'qaccver'.")

    def __iter__(self):
        """Iterates over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()

        if not self._has_comments:
            qresult_key = None
            key_idx = self._key_idx
            while True:
                # get end offset here since we only know a qresult ends after
                # encountering the next one
                end_offset = handle.tell()
                #line = handle.readline()
                line = _bytes_to_string(handle.readline())

                if qresult_key is None:
                    qresult_key = line.split('\t')[key_idx]
                else:
                    try:
                        curr_key = line.split('\t')[key_idx]
                    except IndexError:
                        curr_key = ''

                    if curr_key != qresult_key:
                        yield qresult_key, start_offset, end_offset - start_offset
                        qresult_key = curr_key
                        start_offset = end_offset

                # break if we've reached EOF
                if not line:
                    break
        else:
            # mark of a new query
            query_mark = None
            # mark of the query's ID
            qid_mark = '# Query: '

            while True:
                end_offset = handle.tell()
                line = _bytes_to_string(handle.readline())

                if query_mark is None:
                    query_mark = line
                    start_offset = end_offset
                elif line.startswith(qid_mark):
                    qresult_key = line[len(qid_mark):].split()[0]
                elif line == query_mark or 'BLAST processed' in line:
                    yield qresult_key, start_offset, end_offset - start_offset
                    start_offset = end_offset
                elif not line:
                    break

    def get_raw(self, offset):
        """Returns the raw string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        qresult_raw = ''

        if not self._has_comments:
            key_idx = self._key_idx
            qresult_key = None
            while True:
                line = _bytes_to_string(handle.readline())
                # get the key if the first line (qresult key)
                if qresult_key is None:
                    qresult_key = line.split('\t')[key_idx]
                else:
                    try:
                        curr_key = line.split('\t')[key_idx]
                    except IndexError:
                        curr_key = ''
                    # only break when qresult is finished (key is different)
                    if curr_key != qresult_key:
                        break
                # append to the raw string as long as qresult is the same
                qresult_raw += line
        else:
            # query mark is the line marking a new query
            # something like '# TBLASTN 2.2.25+'
            query_mark = None
            while True:
                line = _bytes_to_string(handle.readline())
                if query_mark is None:
                    query_mark = line
                # if we've encountered another query mark, it's the start of
                # another query
                # if 'BLAST processed' is in line, it's one line before EOF
                elif line == query_mark:
                    break
                # append to the raw string as long as qresult is the same
                qresult_raw += line

                if 'BLAST processed' in line:
                    break

        return _as_bytes(qresult_raw)


class BlastTabWriter(object):

    """Writer for blast-tab output format."""

    def __init__(self, handle, comments=False, fields=_DEFAULT_FIELDS):
        self.handle = handle
        self.has_comments = comments
        self.fields = fields

    def write_file(self, qresults):
        """Writes to the handle, returns how many QueryResult objects are written."""
        handle = self.handle
        qresult_counter, hit_counter, hsp_counter = 0, 0, 0

        for qresult in qresults:
            if self.has_comments:
                handle.write(self.build_comments(qresult))
            if qresult:
                handle.write(self.build_rows(qresult))
                if not self.has_comments:
                    qresult_counter += 1
                hit_counter += len(qresult)
                hsp_counter += sum([len(hit) for hit in qresult])
            # if it's commented and there are no hits in the qresult, we still
            # increment the counter
            if self.has_comments:
                qresult_counter += 1

        # commented files have a line saying how many queries were processed
        if self.has_comments:
            handle.write('# BLAST processed %i queries' % qresult_counter)

        return qresult_counter, hit_counter, hsp_counter

    def build_rows(self, qresult):
        """Returns a string containing tabular rows of the QueryResult object."""
        coordinates = set(['qstart', 'qend', 'sstart', 'send'])
        qresult_lines = ''
        for hit in qresult:
            for hsp in hit:
                line = []
                for field in self.fields:
                    # get the column value ~ could either be an attribute
                    # of qresult, hit, or hsp
                    if field in _COLUMN_QRESULT:
                        value = getattr(qresult, _COLUMN_QRESULT[field][0])
                    elif field in _COLUMN_HIT:
                        value = getattr(hit, _COLUMN_HIT[field][0])
                    # special case, since 'frames' can be determined from
                    # query frame and hit frame
                    elif field == 'frames':
                        value = '%i/%i' % (hsp.query_frame, hsp.hit_frame)
                    elif field in _COLUMN_HSP:
                        value = getattr(hsp, _COLUMN_HSP[field][0])
                    elif field in _COLUMN_FRAG:
                        value = getattr(hsp, _COLUMN_FRAG[field][0])
                    else:
                        assert field not in _SUPPORTED_FIELDS
                        continue

                    # adjust from and to according to strand, if from and to
                    # is included in the output field
                    if field in coordinates:
                        value = self.adjust_coords(field, value, hsp)
                    # adjust output formatting
                    value = self.adjust_output(field, value)

                    line.append(value)

                hsp_line = '\t'.join(line)
                qresult_lines += hsp_line + '\n'

        return qresult_lines

    def adjust_coords(self, field, value, hsp):
        """Adjusts start and end coordinates according to strand."""
        assert field in ('qstart', 'qend', 'sstart', 'send')
        # determine sequence type to operate on based on field's first letter
        seq_type = 'query' if field.startswith('q') else 'hit'

        strand = getattr(hsp, '%s_strand' % seq_type, None)
        if strand is None:
            raise ValueError("Required attribute %r not found." %
                    ('%s_strand' % (seq_type)))
        # switch start <--> end coordinates if strand is -1
        if strand < 0:
            if field.endswith('start'):
                value = getattr(hsp, '%s_end' % seq_type)
            elif field.endswith('end'):
                value = getattr(hsp, '%s_start' % seq_type) + 1
        elif field.endswith('start'):
            # adjust start coordinate for positive strand
            value += 1

        return value

    def adjust_output(self, field, value):
        """Adjusts formatting of the given field and value to mimic native tab output."""

        # evalue formatting, adapted from BLAST+ source:
        # src/objtools/align_format/align_format_util.cpp#L668
        if field == 'evalue':
            if value < 1.0e-180:
                value = '0.0'
            elif value < 1.0e-99:
                value = '%2.0e' % value
            elif value < 0.0009:
                value = '%3.0e' % value
            elif value < 0.1:
                value = '%4.3f' % value
            elif value < 1.0:
                value = '%3.2f' % value
            elif value < 10.0:
                value = '%2.1f' % value
            else:
                value = '%5.0f' % value

        # pident and ppos formatting
        elif field in ('pident', 'ppos'):
            value = '%.2f' % value

        # evalue formatting, adapted from BLAST+ source:
        # src/objtools/align_format/align_format_util.cpp#L723
        elif field == 'bitscore':
            if value > 9999:
                value = '%4.3e' % value
            elif value > 99.9:
                value = '%4.0d' % value
            else:
                value = '%4.1f' % value

        # everything else
        else:
            value = str(value)

        return value

    def build_comments(self, qres):
        """Returns a string of a QueryResult tabular comment."""
        comments = []
        # inverse mapping of the long-short name map, required
        # for writing comments
        inv_field_map = dict((value, key) for key, value in \
                _LONG_SHORT_MAP.items())

        # try to anticipate qress without version
        if not hasattr(qres, 'version'):
            program_line = '# %s' % qres.program.upper()
        else:
            program_line = '# %s %s' % (qres.program.upper(), qres.version)
        comments.append(program_line)
        # description may or may not be present, so we'll do a try here
        try:
            comments.append('# Query: %s %s' % (qres.id, qres.description))
        except AttributeError:
            comments.append('# Query: %s' % qres.id)
        # try appending RID line, if present
        try:
            comments.append('# RID: %s' % qres.rid)
        except AttributeError:
            pass
        comments.append('# Database: %s' % qres.target)
        # qresults without hits don't show the Fields comment
        if qres:
            comments.append('# Fields: %s' % \
                    ', '.join([inv_field_map[field] for field in self.fields]))
        comments.append('# %i hits found' % len(qres))

        return '\n'.join(comments) + '\n'


def _test():
    """Run the Bio.SearchIO.BlastIO.blasttab module's doctests.

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
