#! /usr/bin/env python

formats = ['csv','xml','html']

dr1_url='http://skyserver.sdss.org/dr1/en/tools/search/x_sql.asp'
dr2_url='http://skyserver.sdss.org/dr2/en/tools/search/x_sql.asp'
dr4_url='http://cas.sdss.org/astro/en/tools/search/x_sql.asp'
dr5_url='http://cas.sdss.org/dr5/en/tools/search/x_sql.asp'
dr7_url='http://cas.sdss.org/dr7/en/tools/search/x_sql.asp'
#dr8_url='http://cas.sdss.org/dr8/en/tools/search/x_sql.asp'
dr8_url='http://skyserver.sdss3.org/dr8/en/tools/search/x_sql.asp'
dr9_url='http://skyserver.sdss3.org/dr9/en/tools/search/x_sql.asp'
dr10_url='http://skyserver.sdss3.org/public/en/tools/search/sql.aspx'

default_url=dr9_url
default_fmt='csv'

timeout=30

def usage(status, msg=''):
    "Error message and usage"
    print __doc__
    if msg:
        print '-- ERROR: %s' % msg
    sys.exit(status)

def filtercomment(sql):
    "Get rid of comments starting with --"
    import os
    fsql = ''
    for line in sql.split('\n'):
        fsql += line.split('--')[0] + ' ' + os.linesep;
    return fsql

def query(sql,url=default_url,fmt=default_fmt):
    "Run query and return file object"
    import urllib,socket
    socket.setdefaulttimeout(timeout)
    fsql = filtercomment(sql)
    params = urllib.urlencode({'cmd': fsql, 'format': fmt})
    try:
    	return urllib.urlopen(url+'?%s' % params)    
    except:
    	return urllib.urlopen(url) 

def write_header(ofp,pre,url,qry):
    import  time
    ofp.write('%s SOURCE: %s\n' % (pre,url))
    ofp.write('%s TIME: %s\n' % (pre,time.asctime()))    
    ofp.write('%s QUERY:\n' % pre)
    for l in qry.split('\n'):
        ofp.write('%s   %s\n' % (pre,l))
    
def main(argv):
    "Parse command line and do it..."
    import os, getopt, string
    
    queries = []
    url = os.getenv("SQLCLURL",default_url)
    fmt = default_fmt
    writefirst = 1
    verbose = 0
    
    # Parse command line
    try:
        optlist, args = getopt.getopt(argv[1:],'s:f:q:vlh?')
    except getopt.error, e:
        usage(1,e)
        
    for o,a in optlist:
        if   o=='-s': url = a
        elif o=='-f': fmt = a
        elif o=='-q': queries.append(a)
        elif o=='-l': writefirst = 0
        elif o=='-v': verbose += 1
        else: usage(0)
        
    if fmt not in formats:
        usage(1,'Wrong format!')

    # Enqueue queries in files
    for fname in args:
        try:
            queries.append(open(fname).read())
        except IOError, e:
            usage(1,e)

    # Run all queries sequentially
    for qry in queries:
        ofp = sys.stdout
        if verbose:
            write_header(ofp,'#',url,qry)
        file = query(qry,url,fmt)
        # Output line by line (in case it's big)
        line = file.readline()
        if line.startswith("ERROR"): # SQL Statement Error -&gt; stderr
            ofp = sys.stderr
        if writefirst:
            ofp.write(string.rstrip(line)+os.linesep)
        line = file.readline()
        while line:
            ofp.write(string.rstrip(line)+os.linesep)
            line = file.readline()


if __name__=='__main__':
    import sys
    main(sys.argv)


