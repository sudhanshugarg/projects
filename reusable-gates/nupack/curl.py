# -*- coding: utf-8 -*-

#from lxml import html
from bs4 import BeautifulSoup
import requests
import sys
import re

def main():
#    url = 'http://www.nupack.org/design/new_from_job/61464?token=tz7nWzW292'
#    r = requests.get(url)
    #f = open('curloutput.html','w')
    #f.write(r.text.encode('utf-8'))
    #f.close()
 #   print r.content


    #random = "<html> <body><albums> <rock> <title>Machine Head</title> <artist>Deep Purple</artist> </rock> <blues> <title>Greens From The Garden</title> <artist>Cory Harris</artist> </blues> <country> <title>The Ranch</title> <artist>The Ranch</artist> </country> </albums> </body> </html>"

    #album = BeautifulSoup(random, 'html.parser')
#    album = html.fromstring(random)
#    tree = html.fromstring(r.text.encode('utf-8'))
    #tree = BeautifulSoup(r.text, 'html.parser')
    #tables = tree.xpath('html/body/table/tbody/tr/td/div[@id="main"]/div[@id="design_group"]/table/tbody/tr/td/table/tbody/tr/td/form/table/tbody/tr/td/tt/text()')
    #tables = tree.xpath('//body/table/tbody/tr/td/div[@id="main"]/div[@id="design_group"]/table/tbody/tr/td/table/tbody/tr/td/form/table/tbody/tr/td/tt/text()')
#    tables = tree.xpath('//*[@id="analyze_form_0_0"]/table/tbody/tr[2]/td[1]')

#    titles = album.xpath('html/body/albums/rock/title')
#    print album.albums.rock.title

    #read nupack script from this file.
    with open ('t5d.np', 'r') as content_file:
        design = content_file.read()

    multipart = {
                'design_job[target_structure]': open('t5d.np','rb')
            }
    payload = {
            'preview_token' : '',
            'design_job[nucleic_acid_type]' : 'DNA',
            'design_job[temperature]' : '25.0',
            'design_job[number_of_trials]' : '10',
            'design_job[rna_parameter_file]' : 'rna1995',
            'design_job[dna_parameter_file]' : 'dna1998',
            'design_job[dangle_level]' : '2',
            'design_job[na_salt]' : '0.5',
            'design_job[mg_salt]' : '0.0125',
            'design_job[dotplot_target]' : '0',
            'design_job[prevented_strings]' : 'SSSS\nAAAAA\nAAAA\nGGG',
            'design_job[email_address]' : 'sgarg@cs.duke.edu',
            'commit' : 'Design'
    }
    #print len(design)
    submit_page1 = 'http://www.nupack.org/design/new'
    submit_page2 = 'http://cgi.cs.duke.edu/~sgarg/nupack.php'
    submit_page3 = 'http://localhost:8080'
    hdrs = {
#            'Cache-Control':'max-age=0',
#            'Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
#            'Origin':'null',
#            'User-Agent':'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36',
#            'Content-Type':'multipart/form-data; boundary=----WebKitFormBoundaryRvFIpcbpreh4yDoF',
#            'DNT':'1',
#            'Accept-Encoding':'gzip, deflate',
#            'Accept-Language':'en-US,en;q=0.8',
#            'Cookie':' bs_t_102996039541ecf743394915d27a0cc6=YTo2OntzOjE6InMiO3M6MzI6IjEwMjk5NjAzOTU0MWVjZjc0MzM5NDkxNWQyN2EwY2M2IjtzOjE6ImMiO3M6MzI6ImVkODZhNTg3MTM5OGE5NjgyZjQyOTlmMzQ5OGU2OTFiIjtzOjE6ImQiO3M6MzI6IjhiZDU0MmM1NTM3YjQxMTdhMjMzMWEzODQ4MzQ0N2M5IjtzOjM6ImNpZCI7czoyOToiZHNvdG56a2Q0YjlnMDc5NXl5cjZ0bWhiaW9uYWgiO3M6MzoidGlkIjtzOjY3OiIzLktSQS5DTDZQQlEuRERnZy5BY1QtR0EuLkFwanFhQS5iLi5sLkFST0NrQS5iLlZRZFB6US5WUWVCQlEubGxvamJ3IjtzOjE6ImUiO2k6MTQyNzE0NjM4NDt9;__unam=cc0a4ec-14d06da603d-8421892-4; TWIKISID=b55a28710cb38e9a8f03ce0f82654585'
            }
    rd = requests.post(submit_page1, data=payload, headers=hdrs, files=multipart)
    print rd.text.encode('utf-8')
    #print rd.headers
    #print rd.status_code
'''
    tree = BeautifulSoup(rd.content, 'html.parser')
    tabletag = tree.table.tr.td.find_all(id="main")[0].find_all(id="design_group")[0].table.tr.td.table
    names = tabletag.find_all(string=re.compile('_strand'))
    strands = tabletag.find_all('tt')
    print len(names), names
    print len(strands), strands

    for i in range(len(names)):
        print names[i].strip(), ':', strands[i].string
'''

if __name__ == '__main__':
    main()
