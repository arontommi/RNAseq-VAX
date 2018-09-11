from collections import defaultdict


def read_header(vcfpath):
    """

    get header elements from a vcf file (lines that start with ##)

    """
    with open(vcfpath, 'r') as f:
        headerlines = [l for l in f if l.startswith('##')]
    return headerlines



def clean_header(headerlines):
    """

    cleans headerlines

    """
    cleanheader = []
    for line in headerlines:
        ll = line.replace("\"", "")
        lll = ll.replace(">\n", "")
        llll = lll.replace("\n", "")
        cleanheader.append(llll)
    return cleanheader


def remove_hashtag(string):
    """

    removes ## and returns str

    """
    st = string.replace('##', '')
    return st


def print_header_lines(headerlines):

    """
    creates header structure from headers, returns the header and a sub header if exists
    """

    headers = []
    for line in headerlines:
        if "=<" not in line:
            lineheader = line.split('=')
            headers.append(remove_hashtag(lineheader[0]))
        else:
            cline =line.replace('<', '')
            lineheader = cline.split('=')
            headers.append(tuple([remove_hashtag(lineheader[0]), lineheader[2].split(',')[0]]))
    return headers

def get_last_element(line, spliters):
    """

    remove spliter from line
    """


    for i in spliters:
        if i in line:
            ii = line.split('{}'.format(i))[-1]
        else:
            pass
    try:
        return ii
    except NameError:
        return line



def formatlines(cleanheader, vcftag):
    """

    obtain data for each header

    """
    if isinstance(vcftag, tuple):
        for i in cleanheader:
            if i.startswith('##{0}=<ID={1},'.format(vcftag[0],vcftag[1])):
                return get_last_element(i,  ['length=', 'Description=', 'Format:'])
    else:
        for i in cleanheader:
            if i.startswith('##{0}'.format(vcftag)):
                return get_last_element(i, ['length=','Description=', 'Format:'])


def populate_dict(cleanheader, headers):
    """

    create a dict from header

    """
    d = defaultdict(dict)
    for i in headers:
        if isinstance(i, tuple):
            d[i[0]][i[1]] = formatlines(cleanheader, i)
        else:
            d[i] = formatlines(cleanheader, i)
    return d


def get_header_dict(inputfile):
    """

    main function to obtain headers

    :param inputfile: vcf file
    :return: header dict
    """
    headerlines = read_header(inputfile)
    cleanheader = clean_header(headerlines)
    headers = print_header_lines(cleanheader)
    header_dict = populate_dict(cleanheader, headers)
    return header_dict


def read_csq(CSQ_header_subdict):
    """

    :param CSQ_header_subdict:  uses subdict (dict['INFO']['CSQ']
    :return: list with csq elements

    """
    csq_holder = CSQ_header_subdict.split('|')
    return csq_holder