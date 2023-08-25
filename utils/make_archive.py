#!/usr/bin/env python

"""
Simple script to make .zip archives of this repository, suitable for
upload to Zenodo or a similar DOI-providing service. Unlike Zenodo's
own GitHub integration (or 'git archive') the repository is split into
several smaller archives, so that users don't have to download enormous
archives containing trajectories just to get a small file like a Python
script or alignment.
"""

from __future__ import print_function, division
import sys
import os
import hashlib
import shutil
import subprocess

# Put larger directories (keys) in uniquely-named zipfiles (values)
top_dir = '/home/ignacia/Research/yeast_NPC/modeling_2023/repo/'
ARCHIVES = {
  os.path.join(top_dir,'data'): 'data',
  os.path.join(top_dir,'scripts'): 'scripts',
  os.path.join(top_dir,'results'): 'models'}

REPO="NPC_TMDs"

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def format_size(b):
    suffixes = ['B', 'KiB', 'MiB', 'GiB']
    for s in suffixes:
        if b < 1024:
            return "%.2f %s" % (b, s)
        b /= 1024
    return "%.2f TiB" % b

class Archiver(object):
    ARCHIVE_DIR = "for_archival"

    def __init__(self, tag):
        if os.path.exists(self.ARCHIVE_DIR):
            raise ValueError("The %s directory already exists - please "
                             "delete it first" % self.ARCHIVE_DIR)
        self.tag = tag
        os.mkdir(self.ARCHIVE_DIR)
        self.topdir = os.path.join(self.ARCHIVE_DIR, "%s-%s" % (REPO, self.tag))

    def get_all_files(self):
        print("Extracting all files from %s at %s" % (REPO, self.tag))
        subprocess.check_call('git archive --format=tar --prefix=util/%s/ %s '
                              '| tar -xf -' % (self.topdir, self.tag),
                              shell=True, cwd='..')

    def zip_subdir(self, subdir, zipname):
        base = os.path.basename(subdir)
        print('base', base)
        subdir_full = os.path.join(self.topdir, subdir)
        cwd = os.getcwd()
        print("Archiving %s" % subdir)
        outzip_full = os.path.join(cwd, self.ARCHIVE_DIR, zipname)
        print('aaa', outzip_full, base)
        subprocess.call(['zip', '-r', outzip_full, subdir])
                        #cwd=os.path.join(subdir_full, '..'))
        print('bbb')
        #shutil.rmtree(subdir_full)
        #os.mkdir(subdir_full)
        print('ccc')
        #with open(os.path.join(self.ARCHIVE_DIR, 'README.txt'), 'w') as fh:
        #    fh.write("""
#The files in this directory can be found in the %s file
#at the same DOI where this archive is available.
#""" % zipname)

    def zip_toplevel(self):
        print("Archiving top level")
        dirname = "%s-%s" % (REPO, self.tag)
        subprocess.check_call(['zip', '-r', dirname + '.zip', dirname],
                              cwd=self.ARCHIVE_DIR)
        shutil.rmtree(self.topdir)

    def summarize(self):
        for fname in sorted(os.listdir(self.ARCHIVE_DIR)):
            fullname = os.path.join(self.ARCHIVE_DIR, fname)
            sz = os.stat(fullname).st_size
            print("%s %-10s %s" % (md5(fullname), format_size(sz), fname))
        print("zip files created in %s. Upload them and then"
              % self.ARCHIVE_DIR)
        print("delete that directory.")


def main():
    if len(sys.argv) != 2:
        print("Usage: %s tag" % sys.argv[0], file=sys.stderr)
        sys.exit(1)

    tag = sys.argv[1]
    a = Archiver(tag)
    a.get_all_files()
    # Sort dirs by length, longest first, so we zip any subdirs first
    for subdir in sorted(ARCHIVES.keys(), key=lambda a: -len(a)):
        zipname = ARCHIVES[subdir]
        print('-------',subdir)
        a.zip_subdir(subdir, zipname + '.zip')
    #a.zip_toplevel()
    a.summarize()

if __name__ == '__main__':
    main()
