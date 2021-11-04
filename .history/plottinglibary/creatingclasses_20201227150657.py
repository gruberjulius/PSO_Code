
class apple:
    def printname():
        print(self.i)
    def init(self, i):
        self.i = str(i)



def main():
    numberofclasses = int(input("How many classes do you want to create? "))
    i = 0
    while i  < numberofclasses:
        x='apple' + str(i)   = apple(i)
        exec("%s.apple.printname()" %(x))
    
if __name__ == '__main__':
    main()