
class apple:
    
    def __init__(self, i):
        self.i = i
    
    def printname():
        print(self.i)



def main():
    numberofclasses = int(input("How many classes do you want to create? "))
    i = 0
    while i  < numberofclasses:
        x='apple' + str(i)   #= apple(i)
        x = apple(i)
        #exec("%s.init(%d)"% (x, i))
        #exec("%s.apple.printname()" %(x))
    
    
    
if __name__ == '__main__':
    main()