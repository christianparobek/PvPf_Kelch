'''
Created on Feb 11, 2014
edited on June 15, 2015
@author: ddeconti
edited by Nbrazeau for vennoutput July 20 
'''

from scipy import stats
import getopt
import numpy
import pysam
import sys
import json
import pprint 


def get_major_allele(allele_dict):
    count_dict = {'a': allele_dict['f']['a'] + allele_dict['r']['a'], \
                  't': allele_dict['f']['t'] + allele_dict['r']['t'], \
                  'c': allele_dict['f']['c'] + allele_dict['r']['c'], \
                  'g': allele_dict['f']['g'] + allele_dict['r']['g']}
    return max(count_dict.iterkeys(), key=lambda k: count_dict[k])

def make_vector(allele, allele_dict):
    return [allele_dict['f'][allele], allele_dict['r'][allele]]

def is_strand_bias(allele, allele_dict, min_num=1):
    is_bias = True
    if allele_dict['f'][allele] > min_num and \
    allele_dict['r'][allele] > min_num:
        is_bias = False
    return is_bias

def parse_pileup(bam_filename, ref_filename, min_map_qual=10, \
                 min_read_qual=35, maxdepth=100000, \
                 flag_read_ends=None, flag_dups=False, \
                 near_qual=None, out_filename=None, min_read_size=None, \
                 begin=None, end=None, primer_size=None, \
                 pass_both_strands=True, cumulative_error=False, \
                 skip_indel_reads=False, relative_to_ref=False, aln_score=None):
    total_out_list = []
    try:
        bam_file = pysam.Samfile(bam_filename, 'rb')
        ref_file = pysam.Fastafile(ref_filename)
    except IOError as err:
        sys.exit()
    if out_filename:
        out_venn=out_filename+"venntable.txt"
        venn_handle=open(out_venn, 'w')
        try:
            handle = open(out_filename, 'w')
        except IOError as err:
            sys.stderr.write("IOError: ", str(err) + "\n")
            sys.stderr.flush()
            sys.exit()
    out_list = ["mapq=" + str(min_map_qual), \
                "readq=" + str(min_read_qual), \
                "read_ends=" + str(flag_read_ends), \
                "neighbordis=" + str(near_qual), \
                "min_read_size=" + str(min_read_size), \
                "primer_size=" + str(primer_size), \
                "both_strands=" + str(pass_both_strands), \
                "cumulative_error=" + str(cumulative_error), \
                "indel_f=" + str(skip_indel_reads), \
                "total_area=" + str((end-primer_size) - (begin+primer_size))]
    outstr = '#' + ','.join(out_list) + '\n'
    out_list = ["chrom", \
                "pos", \
                "ref_allele", \
                "major_allele", \
                "minor_allele", \
                "filtered_depth", \
                "minor_allele_freq", \
                "major_f", \
                "major_r", \
                "minor_f", \
                "minor_r", \
                "unfiltered_depth", \
                "mapq_f", \
                "read_size_f", \
                "readq_f", \
                "neighbor_f", \
                "read_ends_f", \
                "indel_f", \
                "bias_f", \
                "as_f", \
                "dup_f", \
                "total_median_depth", \
                "minor_median_depth"]
    outstr += '#' + '\t'.join(out_list) + '\n'
    if out_filename:
        handle.write(outstr)
    else:
        sys.stdout.write(outstr)
        sys.stdout.flush()
    total_depth = []
    minor_allele_depth = []
    for pileup in bam_file.pileup(max_depth=maxdepth):
        if begin != None and end != None and primer_size != None:
            if pileup.pos < begin + primer_size or \
            pileup.pos > end - primer_size:
                    continue
        allele_dict = {'f': {'a': 0, \
                             't': 0, \
                             'c': 0, \
                             'g': 0},\
                       'r': {'a': 0, \
                             't': 0, \
                             'c': 0, \
                             'g': 0}}
        filtered_depth = 0
        full_depth = 0
        filtered_by_mapq = 0
        filtered_by_readq = 0
        filtered_by_read_ends = 0
        filtered_by_dups = 0
        filtered_by_neighbor_readq = 0
        filtered_by_read_size = 0
        filtered_by_indel = 0
        filtered_by_del = 0
        filtered_by_bias = {'a': 0, 't': 0, 'c': 0, 'g': 0}
        filtered_by_as = 0
        ref_allele = ref_file.fetch(bam_file.getrname(pileup.tid),
                                    pileup.pos, pileup.pos+1)
        if flag_dups:
            filtered_pileup_list = []
        
        
        filt_dict={'00000':0, '00001':0, '00010':0, '00011':0, '00100':0, '00101':0, '00110':0, '00111':0,'01000':0,'01001':0,'01010':0,'01011':0,'01100':0,'01101':0,'01110':0,'01111':0,'10000':0,'10001':0,'10010':0,'10011':0,'10100':0,'10101':0,'10110':0,'10111':0,'11000':0,'11001':0, '11010':0, '11011':0,'11100':0,'11101':0,'11110':0, '11111':0}
        
        for pileup_read in pileup.pileups:
            filt_mapq=False
            filt_readq=False
            filt_read_ends=False
            filt_readsize=False
            filt_as=False
        
        
            skip_read = False
            full_depth += 1
            if skip_indel_reads:
                has_indel = False
                for c in pileup_read.alignment.cigar:
                    if c[0] in (1,2):
                        has_indel = True
                        filtered_by_indel += 1
                        break
                if has_indel:
                    if cumulative_error:
                        skip_read = True
                    else:
                        continue
            if pileup_read.is_del:
                filtered_by_del += 1
                if cumulative_error:
                    skip_read = True
                else:
                    continue
            if pileup_read.alignment.mapq < min_map_qual:
                filtered_by_mapq += 1
                filt_mapq=True
                if cumulative_error:
                    skip_read = True
                else:
                    continue
            if pileup_read.alignment.rlen < min_read_size:
                filtered_by_read_size += 1
                filt_readsize=True
                if cumulative_error:
                    skip_read = True
                else:
                    continue
            if aln_score and pileup_read.alignment.tags[0][1] < aln_score:
                filtered_by_as += 1
                filt_as=True
                if cumulative_error:
                    skip_read = True
                else:
                    continue
            if pileup_read.indel == 0:
                base = pileup_read.alignment.seq[pileup_read.query_position].lower()
                qual = ord(pileup_read.alignment.qual[pileup_read.query_position])-33
                if qual < min_read_qual:
                    filtered_by_readq += 1
                    filt_readq=True
                    if cumulative_error:
                        skip_read = True
                    else:
                        continue
                if near_qual:
                    to_skip = False
                    for i in range(1, near_qual+1):
                        p_aln = pileup_read.alignment
                        if p_aln.rlen-i < 0:
                            rq = ord(p_aln.qual[pileup_read.query_position+i])-33
                            if rq < min_read_qual:
                                to_skip = True
                                break
                        elif p_aln.rlen+i > p_aln.rlen-1:
                            lq = ord(p_aln.qual[pileup_read.query_position-i])-33
                            if lq < min_read_qual:
                                to_skip = True
                                break
                        else:
                            lq = ord(p_aln.qual[pileup_read.query_position-i])-33
                            rq = ord(p_aln.qual[pileup_read.query_position+i])-33
                            if lq < min_read_qual or rq < min_read_qual:
                                to_skip = True
                                break
                    if to_skip:
                        filtered_by_neighbor_readq += 1
                        if cumulative_error:
                            skip_read = True
                        else:
                            continue
                if flag_read_ends:
                    if pileup_read.alignment.rlen - pileup_read.query_position < \
                    flag_read_ends:
                        filtered_by_read_ends += 1
                        filt_read_ends=True
                        if cumulative_error:
                            skip_read = True
                        else:
                            continue
                    if pileup_read.query_position < flag_read_ends:
                        filtered_by_read_ends += 1
                        if cumulative_error:
                            skip_read = True
                        else:
                            continue
                if cumulative_error and skip_read:
                    bstr="{0:b}{1:b}{2:b}{3:b}{4:b}".format(filt_mapq, filt_readq, filt_read_ends, filt_readsize, filt_as)
                    my_list=''.join(bstr)
                    if my_list in filt_dict:
                        filt_dict[my_list]+=1
                    else:
                        sys.stderr.write("Combination not in dictionary")
                        sys.exit(1)
                else:
                    continue
                if pileup_read.alignment.is_reverse:
                    allele_dict['r'][base] += 1
                else:
                    allele_dict['f'][base] += 1
                if flag_dups:
                    filtered_pileup_list.append(pileup_read)
                filtered_depth += 1
        major_allele = get_major_allele(allele_dict)
        
        for key, value in filt_dict.items():
            venn_handle.write(str([key, value])+'\n')
        
        # duplicate removal
        # Experimental
        if flag_dups:
            minor_allele_keys = ['a', 't', 'c', 'g']
            minor_allele_keys.remove(major_allele)
            minor_dict = {}
            for minor_allele in minor_allele_keys:
                minor_dict[minor_allele] = {'f': {}, 'r': {}}
               # minor_dict[minor_allele] = {'f': None, 'r': None}
            for pileup_read in filtered_pileup_list:
                base = pileup_read.alignment.seq[pileup_read.query_position].lower()
                if base not in minor_allele_keys:
                    continue
                read_pos_hash = (pileup_read.alignment.pos, \
                             pileup_read.alignment.aend)
                if pileup_read.alignment.is_reverse:
                    dir_key = 'r'
                else:
                    dir_key = 'f' 
                #if not minor_dict[base][dir_key]:
                #    minor_dict[base][dir_key] = {read_pos_hash: [pileup_read]}
                if read_pos_hash in minor_dict[base][dir_key]:
                    minor_dict[base][dir_key][read_pos_hash].append(pileup_read)
                else:
                    minor_dict[base][dir_key][read_pos_hash] = [pileup_read]
            #print "before:", allele_dict
            filtered_depth = 0
            for base in minor_allele_keys:
                if len(minor_dict[base]['f'].keys()) < flag_dups and \
                len(minor_dict[base]['r'].keys()) < flag_dups:
                    allele_dict['f'][base] = 0
                    allele_dict['r'][base] = 0
                    #allele_dict['f'][base] = len(minor_dict[base]['f'].keys())
                    #filtered_depth += len(minor_dict[base]['f'].keys()) + \
                    #len(minor_dict[base]['r'].keys())
                    #allele_dict['r'][base] = len(minor_dict[base]['r'].keys())
                else:
                    filtered_depth += allele_dict['f'][base] + \
                    allele_dict['r'][base]
            filtered_depth += allele_dict['f'][major_allele] + \
            allele_dict['r'][major_allele]
            #print "after:", allele_dict
        # end duplicate removal section
        
        total_depth.append(filtered_depth)
        passed_first = False
        #major_vector = make_vector(major_allele, allele_dict)
        for base in ('a', 't', 'c', 'g'):
            is_pass = 0
            if base == major_allele:
                if not relative_to_ref or base == ref_allele.lower():
                    continue
            if pass_both_strands and is_strand_bias(base, allele_dict):
                continue
            if not passed_first:
                minor_allele_depth.append(filtered_depth)
                passed_first = True
            minor_vector = make_vector(base, allele_dict)
            temp_vector_list = []
            for b in ('a', 't', 'c', 'g'):
                if b == base:
                    continue
                else:
                    temp_vector_list.append(make_vector(b, allele_dict))
            other_vector = [0,0]
            for temp_vector in temp_vector_list:
                other_vector = numpy.add(other_vector, temp_vector)
            base_count = numpy.sum(minor_vector)
            try:
                minor_allele_freq = base_count/float(filtered_depth)
            except ZeroDivisionError:
                minor_allele_freq = 0
            #oddsratio, pvalue = stats.fisher_exact([major_vector, \
            #                                        minor_vector])
            out_list = [bam_file.getrname(pileup.tid), \
                        str(pileup.pos), \
                        ref_allele, \
                        major_allele, \
                        base, \
                        str(filtered_depth), \
                        str(minor_allele_freq), \
                        str(other_vector[0]), \
                        str(other_vector[1]), \
                        #str(allele_dict['f'][major_allele]), \
                        #str(allele_dict['r'][major_allele]), \
                        str(allele_dict['f'][base]), \
                        str(allele_dict['r'][base]), \
                        str(full_depth), \
                        str(filtered_by_mapq), \
                        str(filtered_by_read_size), \
                        str(filtered_by_readq), \
                        str(filtered_by_neighbor_readq), \
                        str(filtered_by_read_ends), \
                        str(filtered_by_indel), \
                        str(filtered_by_bias[base]), \
                        str(filtered_by_as), \
                        str(filtered_by_dups)]
            total_out_list.append(out_list)
            #if out_filename:
            #    handle.write('\t'.join(out_list) + '\n')
            #    #sys.stdout.write('\t'.join(out_list) + '\n')
            #    #sys.stdout.flush()
            #else:
            #    sys.stdout.write('\t'.join(out_list) + '\n')
            #    sys.stdout.flush()
    if len(total_depth) > 0:
        total_median = int(round(numpy.median(total_depth), 0))
    else:
        total_median = 0
    if len(minor_allele_depth) > 0:
        minor_allele_median = int(round(numpy.median(minor_allele_depth), 0))
    else:
        minor_allele_median = 0
    for out_list in total_out_list:        
        if out_filename:
            out_list.append(str(total_median))
            out_list.append(str(minor_allele_median))
            handle.write('\t'.join(out_list) + '\n')
        else:
            out_list.append(str(total_median))
            out_list.append(str(minor_allele_median))
            sys.stdout.write('\t'.join(out_list) + '\n')
            sys.stdout.flush()

def usage():
    outstr = "\nUsage: minor_allele_catcher [options]" + \
    "<ref.fasta> <in.sorted.bam>\n\n" + \
    "If no default given, either false or 0\n" + \
    "Options: -q/--qual           minimum read qual [default: 10]\n" + \
    "         -m/--mapq           minimum map qual [default: 10]\n" + \
    "         -e/--ends           minimum distance from read ends\n" + \
    "         -d/--depth          max depth for pysam.pileup" + \
    " [default: 100000\n" + \
    "         -r/--rmdup          remove duplicates\n" + \
    "         -n/--nearq          minimum distance near snp of passing" + \
    " quality\n" +\
    "         -s/--size           minimum read size\n" + \
    "         -b/--bothstrands    filter if strand bias [default: False]\n" +\
    "         -c/--cumulative     count filter cumulatively" +\
    " [default: False]\n" +\
    "         -i/--indel          skip indel reads [default: False]\n" + \
    "         --refrelative       compare snps to ref rather" +\
    "         -a/--aln_score      minimum alignment score\n" + \
    " than major allele [default: false]\n" + \
    "         -p/--primersize     size of primer\n" + \
    "         -l/--lower          minimum alignment range\n" + \
    "         -u/--upper          maximum alignment range\n" + \
    "         -o/--out            output filename\n" + \
    "         -h/--help           help screen\n\n"
    sys.stderr.write(outstr)
    sys.stderr.flush()
    sys.exit()

def main(sa):
    if len(sa) < 2:
        usage()
    bam_filename = sa[-1]
    ref_filename = sa[-2]
    min_map_qual = 10
    min_read_qual = 10
    max_depth = 100000
    flag_read_ends = None
    flag_dups = False
    near_qual = None
    out_filename = None
    min_read_size = None
    begin = 0
    end = 10000000
    primer_size = 0
    pass_both_strands = False
    cumulative_error = False
    skip_indel_reads = False
    relative_to_ref = False
    aln_score = True
    try:
        opt_list, arg_list = getopt.getopt(sa[:-2], \
                                           "q:m:e:d:r:n:o:s:u:l:p:bcia:h", \
                                           ["qual=",\
                                            "mapq=", \
                                            "ends=", \
                                            "depth=", \
                                            "rmdup=", \
                                            "nearq=", \
                                            "out=", \
                                            "size=", \
                                            "upper=", \
                                            "lower=", \
                                            "primersize=", \
                                            "bothstrands", \
                                            "cumulative", \
                                            "indel", \
                                            "aln_score", \
                                            "refrelative", \
                                            "help"])
    except getopt.GetoptError as err:
        sys.stderr.write("GetoptError: " + str(err) + "\n")
        sys.stderr.flush()
        usage()
    for opts, args in opt_list:
        if opts in ("-q", "--qual"):
            try:
                min_read_qual = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("-m", "--mapq"):
            try:
                min_map_qual = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("-e", "--ends"):
            try:
                flag_read_ends = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("-d", "--depth"):
            try:
                max_depth = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("-r", "--rmdup"):
            try:
                flag_dups = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("-n", "--nearq"):
            try:
                near_qual = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("-o", "--out"):
            out_filename = args
        elif opts in ("-s", "--size"):
            try:
                min_read_size = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("-u", "--upper"):
            try:
                end = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("-l", "--lower"):
            try:
                begin = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("-p", "--primersize"):
            try:
                primer_size = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("-b", "--bothstrands"):
            pass_both_strands = True
        elif opts in ("-c", "--cumulative"):
            cumulative_error = True
        elif opts in ("-i", "--indel"):
            skip_indel_reads = True
        elif opts in ("-a", "--aln_score"):
            try:
                aln_score = int(args)
            except ValueError as err:
                sys.stderr.write("ValueError: " + str(err) + "\n")
                usage()
        elif opts in ("--refrelative"):
            relative_to_ref = True
        elif opts in ("-h", "--help"):
            usage()
        else:
            usage()
    if begin > end:
        print begin, end
        sys.stderr.write("Left edge > right edge\n")
        usage()
    parse_pileup(bam_filename, ref_filename, min_map_qual, min_read_qual, \
                 max_depth, flag_read_ends, flag_dups, near_qual, \
                 out_filename, min_read_size, begin, end, primer_size, \
                 pass_both_strands, cumulative_error, skip_indel_reads, \
                 relative_to_ref, aln_score)
    sys.exit()

if __name__ == '__main__':
    main(sys.argv[1:])