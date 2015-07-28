# -*- coding: utf-8 -*-

#from lxml import html
from bs4 import BeautifulSoup
import requests
import sys
import re

def main():
    r = requests.get(sys.argv[1])
    f = open('curloutput.html','w')
    f.write(r.text.encode('utf-8'))
    f.close()


    random = "<html> <body><albums> <rock> <title>Machine Head</title> <artist>Deep Purple</artist> </rock> <blues> <title>Greens From The Garden</title> <artist>Cory Harris</artist> </blues> <country> <title>The Ranch</title> <artist>The Ranch</artist> </country> </albums> </body> </html>"

    album = BeautifulSoup(random, 'html.parser')
#    album = html.fromstring(random)
#    tree = html.fromstring(r.text.encode('utf-8'))
    tree = BeautifulSoup(r.text, 'html.parser')
    #tables = tree.xpath('html/body/table/tbody/tr/td/div[@id="main"]/div[@id="design_group"]/table/tbody/tr/td/table/tbody/tr/td/form/table/tbody/tr/td/tt/text()')
    #tables = tree.xpath('//body/table/tbody/tr/td/div[@id="main"]/div[@id="design_group"]/table/tbody/tr/td/table/tbody/tr/td/form/table/tbody/tr/td/tt/text()')
#    tables = tree.xpath('//*[@id="analyze_form_0_0"]/table/tbody/tr[2]/td[1]')

#    titles = album.xpath('html/body/albums/rock/title')
#    print album.albums.rock.title
    tabletag = tree.table.tr.td.find_all(id="main")[0].find_all(id="design_group")[0].table.tr.td.table
    names = tabletag.find_all(string=re.compile('_strand'))
    strands = tabletag.find_all('tt')
    print len(names), names
    print len(strands), strands

    for i in range(len(names)):
        print names[i].strip(), ':', strands[i].string


if __name__ == '__main__':
    main()
