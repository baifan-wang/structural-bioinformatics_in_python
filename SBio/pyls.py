import os

class ls():
    def __init__(self, directory="."):
        self.dir = directory
        if self.dir is not ".":
            if not os.path.exists(self.dir):
                raise FileNotFoundError("No such file or directory: '%s'" %self.dir)
            if not os.path.isdir(self.dir):
                raise TypeError("'%s' is not a directory" %self.dir)

    def all_files(self):
        return sorted([os.path.join(self.dir,x) for x in os.listdir(self.dir) \
        if os.path.isfile(os.path.join(self.dir,x))])

    def all_specified_files(self, file_extenstion):
        fe = file_extenstion
        return sorted([os.path.join(self.dir,x) for x in os.listdir(self.dir) \
        if os.path.isfile(os.path.join(self.dir,x)) and os.path.splitext(x)[1]==fe])

    def all_folders(self):
        return sorted([os.path.join(self.dir,x) for x in os.listdir(self.dir) \
        if os.path.isdir(os.path.join(self.dir,x))])

    def __getattr__(self, args):
        if args == 'all':
            return self.all_files()
        elif args == 'folder':
            return self.all_folders()
        else:
            fe="."+args              #file extension
            return self.all_specified_files(fe)


if __name__ == "__main__":
    print(ls().folder)
    print(ls().all)