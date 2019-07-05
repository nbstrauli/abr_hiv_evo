import os

class master_directory(object):
    """
    A class working with directories that have an ambiguous file structure within them

    Attributes:
    dirpath - path to the master directory in question.
    file_count - Number of files (non-directories) in the directory tree. These do not include hidden files, such as those that start with '.'.
    dirs_containing_filetype - This a set datatype where each element is a path to a directory that contains a provided filetype. This attribute is initially set to be an empty set. Must run 'get_directories_containing_filetype' to populate this attribute.
    files_containing_filetype - A set of absolute filepaths that are of files of a given query type (i.e. fasta, fastq, txt, etc.). This is and empty set until the 'get_filepaths_containing_filetype', or 'get_filepaths_containing_filetypes', or 'mirror_directory_tree_with_files' methods are run.
    mirrored_directory_leaves - This is a set of the paths to each of the directories that contain leaves in the mirrored directory. One must run 'mirror_directory_tree' method for this attribute to be a non-empty set.
    mirrored_filepaths - These are the paths to the analogous files that would populate the mirrored directory tree. They are analogous to the files in the 'files_containing_filetype' attribute. This is an empty set until 'mirror_directory_tree_with_files' is run.
    """
    def __init__(self, dirpath):
        if dirpath[-1] != '/':
            dirpath += '/'
        if not os.path.isdir(dirpath):
            print dirpath, 'either is not a directory or does not exist.'
        self.dirpath = dirpath
        self.file_count = None
        self.dirs_containing_filetype = set()
        self.files_containing_filetype = set()
        self.mirrored_directory_leaves = set()
        self.mirrored_filepaths = set()
        return

    def count_files_loop(self, dirpath):
        """A self referential method. We could do this with os.walk probably a lot easier, but we wanted to try a cool application of a self referential funciton. This method cycles through all the levels of a directory tree and counts the number of files it encounters along the way. Updates that file count in self.file_count."""
        for i in os.listdir(dirpath):
            if i[0] == '.':
                continue
            elif os.path.isdir(dirpath + i):
                self.count_files_loop(dirpath + i + '/')
            elif os.path.isfile(dirpath + i):
                self.file_count += 1
            else:
                print dirpath + i, 'does not exist'
        return

    def count_files(self):
        """Simply runs the 'count_files_loop' above."""
        self.file_count = 0
        self.count_files_loop(self.dirpath)
        return

    @staticmethod
    def get_directories_containing_filetype_loop(self, dirpath, filetype, len_suffix):
        """This method cycles through the entire file structure and retrieves all the directories that contain files with a given suffix (ex: file names ending in '.fastq'). It then returns a set of the paths to these directories."""
        for i in os.listdir(dirpath):
            if i[0] == '.':
                continue
            elif os.path.isdir(dirpath + i):
                self.get_directories_containing_filetype_loop(self, dirpath + i + '/', filetype, len_suffix)
            elif os.path.isfile(dirpath + i):
                if i[-len_suffix:] == filetype:
                    self.dirs_containing_filetype.update([dirpath])
            else:
                print dirpath + i, 'does not exist'
        return
    
    def get_directories_containing_filetype(self, filetype):
        """This method simply runs 'get_directories_containing_filetype_loop' above. However, one should probably use this method as opposed to that above, as the input is simpler."""
        len_suffix = len(filetype)
        self.get_directories_containing_filetype_loop(self, self.dirpath, filetype, len_suffix)
        return

    @staticmethod
    def get_filepaths_containing_filetype_loop(self, dirpath, filetype, len_suffix):
        """This method does the same as 'get_directories_containing_filetype_loop' but gets the filepaths that have a query suffix, rather than the containing directory paths."""
        for i in os.listdir(dirpath):
            if i[0] == '.':
                continue
            elif os.path.isdir(dirpath + i):
                self.get_filepaths_containing_filetype_loop(self, dirpath + i + '/', filetype, len_suffix)
            elif os.path.isfile(dirpath + i):
                if i[-len_suffix:] == filetype:
                    self.files_containing_filetype.update([dirpath + i])
            else:
                print dirpath + i, 'does not exist'
        return

    def get_filepaths_containing_filetype(self, filetype):
        """This method runs the above static method 'get_filepaths_containing_filetype_loop'. It's a bit cleaner to start the self referential loop using this method."""
        len_suffix = len(filetype)
        self.get_filepaths_containing_filetype_loop(self, self.dirpath, filetype, len_suffix)
        return

    @staticmethod
    def get_filepaths_containing_filetypes_loop(self, dirpath, filetypes, len_suffixs, avoid_dirpath=None, avoid_files_with=None):
        """Same as 'get_filepaths_containing_filetype_loop', but will return the filepaths for multiple query filetypes."""
        for i in os.listdir(dirpath):
            if i[0] == '.':
                continue
            elif os.path.isdir(dirpath + i):
                if avoid_dirpath and avoid_dirpath == dirpath+i+'/':
                    continue
                self.get_filepaths_containing_filetypes_loop(self, dirpath + i + '/', filetypes, len_suffixs, avoid_dirpath, avoid_files_with)
            #if is a file and filepath appears to have a suffix
            elif os.path.isfile(dirpath + i) and '.' in dirpath + i:
                if avoid_files_with:
                    file_base = '.'.join(i.split('.')[:-1])
                    if avoid_files_with in file_base:
                        continue
                suffix = i.split('.')[-1]
                if suffix in filetypes:
                    self.files_containing_filetype.update([dirpath + i])
            else:
                print dirpath + i, 'does not exist'
        return

    def get_filepaths_containing_filetypes(self, filetypes, avoid_dirpath=None, avoid_files_with=None):
        """Same as above but can input multiple filetypes. In this case 'filetypes' is a list of strings."""
        len_suffixs = [len(i) for i in filetypes]
        if avoid_dirpath:
            if avoid_dirpath[-1] != '/':
                avoid_dirpath += '/'
            if avoid_dirpath == self.dirpath:
                return
        self.get_filepaths_containing_filetypes_loop(self, self.dirpath, filetypes, len_suffixs, avoid_dirpath, avoid_files_with)
        return

    @staticmethod
    def mirror_directory_tree_loop(self, in_dirpath, out_dirpath, avoid_dirpath=None):
        """This method takes the directory tree in self.dirpath and mirros it in out_dirpath. It does not copy files, but just makes the same directory tree with the same names for the analogous directories. It also updates the self.mirrored_directory_leaves attribute so that it contains a list of the directories that are at all the leaves of the mirrored directory tree. That is, it lists the analogous directories in the mirred tree that contained files (i.e. non-directories) in the self.dirpath directory tree. This ignores README files."""
        for i in os.listdir(in_dirpath):
            if i[0] == '.' or i[:6] == 'README':
                continue
            elif os.path.isdir(in_dirpath + i):
                if avoid_dirpath and avoid_dirpath == in_dirpath+i+'/':
                    continue
                if not os.path.exists(out_dirpath + i):
                    os.makedirs(out_dirpath + i)
                self.mirror_directory_tree_loop(self, in_dirpath + i + '/', out_dirpath + i + '/')
            elif os.path.isfile(in_dirpath + i):
                self.mirrored_directory_leaves.update([out_dirpath])
            else:
                print dirpath + i, 'does not exist'
        return

    def mirror_directory_tree(self, out_dirpath, avoid_dirpath=None):
        """This method runs the 'mirror_directory_tree_loop' above."""
        if out_dirpath[-1] != '/':
            out_dirpath += '/'
        if not os.path.exists(out_dirpath):
            os.makedirs(out_dirpath)
        if avoid_dirpath:
            if avoid_dirpath[-1] != '/':
                avoid_dirpath += '/'
            if self.dirpath == avoid_dirpath:
                return
        self.mirror_directory_tree_loop(self, self.dirpath, out_dirpath, avoid_dirpath)
        return

    @staticmethod
    def mirror_directory_tree_with_files_loop(self, in_dirpath, out_dirpath, only_include_filetypes, include_file_suffix, avoid_files_with):
        """
        This excecutes the method 'mirror_directory_tree_with_files'.
        """
        for i in os.listdir(in_dirpath):
            if i[0] == '.' or i[:6] == 'README':
                continue
            elif os.path.isdir(in_dirpath + i):
                if not os.path.exists(out_dirpath + i):
                    os.makedirs(out_dirpath + i)
                self.mirror_directory_tree_with_files_loop(self, in_dirpath + i + '/', out_dirpath + i + '/', only_include_filetypes, include_file_suffix, avoid_files_with)
            elif os.path.isfile(in_dirpath + i):
                if avoid_files_with:
                    if avoid_files_with in '.'.join(i.split('.')[:-1]):
                        continue
                if only_include_filetypes:
                    suffix = i.split('.')[-1]
                    if suffix in only_include_filetypes:
                        if include_file_suffix:
                            filename = i
                        else:
                            filename = '.'.join(i.split('.')[:-1])
                        self.files_containing_filetype.update([in_dirpath + i])
                        self.mirrored_filepaths.update([out_dirpath + filename])
                        self.mirrored_directory_leaves.update([out_dirpath])
                else:
                    if include_file_suffix or not '.' in i:
                        filename = i
                    else:
                        filename = '.'.join(i.split('.')[:-1])
                    self.files_containing_filetype.update([in_dirpath + i])
                    self.mirrored_filepaths.update([out_dirpath + filename])
                    self.mirrored_directory_leaves.update([out_dirpath])
            else:
                print dirpath + i, 'does not exist'
        return

    def mirror_directory_tree_with_files(self, out_dirpath, only_include_filetypes=None, include_file_suffix=True, avoid_files_with=None):
        """
        This method not only mirrors the directory tree, but it also produces a list of absolute filepaths in the new (mirrored) directory tree that are present in the old directory tree. It does not make new files, it just produces analogous filepaths to the files that are in the original directory.
        out_dirpath - The root of the new (mirrored) directory. Will be at the analogous level to self.dirpath.
        only_include_filetypes - This should be a list of file suffixs (in the form of strings). This will reduce the method to only giving mirrored filepaths for files that have the provided suffixs. For example, if one only wants to mirror fasta formatted files (that end in '.fasta') then this parameter should equal ['fasta']. This way any files in the input directory tree that don't end in 'fasta', will not have a mirrored filepath. One should not include the '.' in the suffixs, just the suffix itself (i.e. should be 'fasta' NOT '.fasta'). If this equals None, then all files will generate an analogous mirrored filepath, regardless of suffix (if any)
        include_file_suffix - If True (default), then the mirrored filepaths will have the suffix of the original filepaths. If False, then the suffix is stripped from the mirrored filepaths.
        avoid_files_with - This means that files whose basename contains this character will no be included.
        """
        if out_dirpath[-1] != '/':
            out_dirpath += '/'
        if not os.path.exists(out_dirpath):
            os.makedirs(out_dirpath)
        self.mirror_directory_tree_with_files_loop(self, self.dirpath, out_dirpath, only_include_filetypes, include_file_suffix, avoid_files_with)
        return
