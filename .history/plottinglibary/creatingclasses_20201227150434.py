
class apple:
    def printname():
        print(self.i)




def main():
    numberofclasses = int(input("How many classes do you want to create? "))
    i = 0
    while i  < numberofclasses:
        x='apple' + str(i)    
        exec("%s.apple.printname()" %(x))
    
if __name__ == '__main__':
    main()