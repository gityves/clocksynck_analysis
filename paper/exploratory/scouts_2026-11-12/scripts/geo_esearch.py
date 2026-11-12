import sys,json,urllib.request,urllib.parse,time
def es(term,n=40):
    u='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmode=json&retmax=%d&term=%s'%(n,urllib.parse.quote(term))
    ids=json.load(urllib.request.urlopen(u))['esearchresult']['idlist']; time.sleep(.4)
    if not ids: return []
    u='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&retmode=json&id='+','.join(ids)
    r=json.load(urllib.request.urlopen(u))['result']; time.sleep(.4)
    return [(r[i]['accession'],r[i]['entrytype'],r[i]['n_samples'],r[i]['title'][:110],r[i].get('gdstype','')[:30]) for i in ids if r[i]['accession'].startswith('GSE')]
for t in sys.argv[1:]:
    print('##',t)
    for x in es(t): print(*x,sep=' | ')
