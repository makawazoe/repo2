#!/usr/bin/python3

import re
import requests
from bs4 import BeautifulSoup

organism_name =  'Lilium longiflorum'

# obtaining taxid
url_1 = 'https://www.ncbi.nlm.nih.gov/taxonomy/?term=' + organism_name

response_1 = requests.get(url_1)
soup_1 = BeautifulSoup(response_1.text, "html.parser")

elems_1 = soup_1.find_all(ref=re.compile("ncbi_uid="))
elems_2 = elems_1[0].attrs['ref']
p3 = r'ncbi_uid=(\d+)'
taxid = re.findall(p3, elems_2)

# obtaining organism info
url_2 = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=' + taxid[0]

response_2 = requests.get(url_2)
soup_2 = BeautifulSoup(response_2.text, "html.parser")

elems_3 = soup_2.body.text
print(elems_3)

#print(elems_1[0].attrs['ref'])
#print(elems_1[0].contents[0].string)
#print(taxid)

#with open("test_4.txt", "a") as f:
#    f.write(str(elems_1))
#    f.write(str(taxid))

with open("test_5.txt", "a") as f:
    f.write(elems_3)
