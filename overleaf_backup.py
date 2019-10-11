
'''
Backup all overleaf projects locally.
Use a git clone to the backup and then loop through pulling.

Adapted from https://github.com/tdurieux/overleaf-backup/blob/master/overleaf_backup.py

To run this, you need:
(1) Change the output script.
(2) Save your email and password for overleaf with git-credentials (https://stackoverflow.com/questions/35942754/how-to-save-username-and-password-in-git).
(3) Clone one project from overleaf to save the credentials (FYI, this will save your password in plain-text in .git-credentials)
(4) Setup a cronjob with `crontab -e` to run the backup every so often.

'''

'''
MIT License

Copyright (c) 2018 Thomas Durieux

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import argparse
import requests
import re
import json
import os
import zipfile
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


osjoin = os.path.join

rootUrl = "https://www.overleaf.com"
pathBackup = os.path.expanduser("~/ownCloud/overleaf_backup/")
isDownloadPDF = False
isDownloadZIP = False


def initParser():
    parser = argparse.ArgumentParser(description='Create a backup of all your overleaf paper')
    parser.add_argument('-email', required=True, help='overleaf email')
    parser.add_argument('-password', required=True, help='overleaf password')
    parser.add_argument('-url', required=False, default=rootUrl, help='The overleaf server')
    parser.add_argument('-output', default=pathBackup, required=False, help='The output folder of the backup')
    parser.add_argument('-pdf', action="store_true", default=False, help='Download the PDF')
    parser.add_argument('-zip', action="store_true", default=False, help='Download the zip')
    return parser.parse_args()


session = requests.Session()


def login(email, password):
    url = rootUrl + "/login"
    r = session.get(url, verify=False)
    m = re.search('window.csrfToken = "([^"]+)";', r.content.decode('utf-8'))
    if m:
        csrf = m.group(1)
        data = {
            "email": email,
            "password": password,
            "_csrf": m.group(1)
        }
        r = session.post(url, data, verify=False)
        return r.status_code == 200
    return False


def getPapers():
    url = rootUrl + "/project"
    r = session.get(url, verify=False)
    m = re.search("\\{\"projects\":([^\n]*)", r.text)
    if m:
        data = json.loads(m.group(0))
        papers = data['projects'][:]
        # for paper in data['projects']:
        #     if paper['archived'] or paper['isV1Project']:
        #         papers.remove(paper)
        return papers
    return []


def updatePaper(papers, raise_on_fail=False):

    git_url = 'https://git.overleaf.com/'

    for paper in papers:

        paper_folder_name = u'' + paper['name'].replace("/", "-").replace(" ", "_").replace(".", "").replace(":", "")

        paper_folder = osjoin(pathBackup, paper_folder_name)

        if not os.path.exists(paper_folder):
            # Clone new projects if not found
            out = os.system(f"git clone {git_url}/{paper['id']} {paper_folder}")
            if out != 0:
                if raise_on_fail:
                    raise ValueError(f"New clone failed for {paper_folder_name}")
                else:
                    print(f"New clone failed for {paper_folder_name}")
            # Don't need to update in this case.
            continue

        curr_dir = os.getcwd()

        # If it exists, force a pull
        os.chdir(paper_folder)
        out = os.system(f"git pull")
        if out != 0:
            # Catch case where this fails.
            # This is a backup that is NOT meant for local changes.
            # So allow a forced reset
            out_force = os.system("git reset --hard HEAD")
            out = os.system(f"git pull")

            if out_force != 0 or out != 0:
                raise ValueError(f"Force pull did not work for {paper_folder_name}")

        os.chdir(curr_dir)


def compile(paper):
    r = session.get(rootUrl + "/project/" + paper['id'], verify=False)
    m = re.search('window.csrfToken = "([^"]+)";', r.content.decode('utf-8'))
    if m:
        csrf = m.group(1)

        url = rootUrl + "/project/" + paper['id'] + "/compile"
        data = {
            "rootDoc_id": paper['id'],
            "_csrf": csrf
        }
        r = session.post(url, data=data, verify=False, stream=True)
        return r.status_code == 200
    return False


def downloadPDF(papers):

    # Make a separate PDFs folder

    pdf_folder = osjoin(pathBackup, 'pdfs')
    if not os.path.exists(pdf_folder):
        os.mkdir(pdf_folder)

    for paper in papers:

        name = u'' + paper['name'].replace("/", "-").replace(" ", "_").replace(".", "").replace(":", "")
        filename = osjoin(pdf_folder, f"{name}_{paper['id']}.pdf".strip())

        if not compile(paper):
            print(f"Failed to compile {paper['name']}.")
            continue

        url = rootUrl + "/project/" + paper['id'] + "/output/output.pdf"
        r = session.get(url, verify=False, stream=True)
        if r.status_code != 200:
            print(f"Failed to connect to {paper['name']}")
            continue

        # Delete the old version if it exists
        if os.path.exists(filename):
            os.remove(filename)

        with open(filename, 'wb') as fd:
            for chunk in r.iter_content(100):
                fd.write(chunk)


def downloadZip(papers):

    # Make a separate PDFs folder

    zip_folder = osjoin(pathBackup, 'zips')
    if not os.path.exists(zip_folder):
        os.mkdir(zip_folder)

    for paper in papers:

        url = rootUrl + "/project/" + paper['id'] + "/download/zip"
        r = session.get(url, verify=False, stream=True)
        if r.status_code != 200:
            print(f"Unable to download zip from {paper['name']}")
            continue

        name = u'' + paper['name'].replace("/", "-").replace(" ", "_").replace(".", "").replace(":", "")

        filename = osjoin(zip_folder, f"{name}_{paper['id']}.zip".strip())

        # Delete the old version if it exists
        if os.path.exists(filename):
            os.remove(filename)

        with open(filename, 'wb') as fd:
            for chunk in r.iter_content(100):
                fd.write(chunk)
        # with open(filename, 'rb') as fh:
        #     z = zipfile.ZipFile(fh)
        #     foldername = os.path.join(pathBackup, paper['id'].strip())
        #     if not os.path.exists(foldername):
        #         os.makedirs(foldername)
        #     for name in z.namelist():
        #         z.extract(name, foldername)
        # os.remove(filename)


if __name__ == '__main__':
    args = initParser()
    rootUrl = args.url
    pathBackup = args.output
    isDownloadPDF = args.pdf
    isDownloadZIP = args.zip

    if login(args.email, args.password):
        papers = getPapers()

        updatePaper(papers)

        if isDownloadPDF:
            downloadPDF(papers)

        if isDownloadZIP:
            downloadZip(papers)

    else:
        print("Invalid password")
