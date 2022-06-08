import os
import sys
import bs4

filepath = sys.argv[1]

with open(filepath, 'rb') as stream:
    soup = bs4.BeautifulSoup(stream.read().decode('utf-8'), 'lxml')
cell_data = []

    # Search out for widgets and place placeholder tag
#for widget in soup.findAll('div', {'class': ['cell_output docutils container']}):
for widget in soup.findAll('script', type='application/vnd.jupyter.widget-state+json'):
    placeholder = ''
    widget.insert_after(placeholder)

stream.close()

html=soup.prettify()
with open(filepath, "wb") as file:
    file.write(html.encode('utf-8'))


###
# 777            <div class="cell_output docutils container">
# 778             <script type="application/vnd.jupyter.widget-view+json">
