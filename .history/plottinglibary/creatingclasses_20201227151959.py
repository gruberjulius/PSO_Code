
class apple:
    
    def __init__(self, i):
        self.i = i
    
    def printname():
        print(self.i)



def main():
    numberofclasses = int(input("How many classes do you want to create? "))
    """
    i = 0
    while i  < numberofclasses:
        x='apple' + str(i)   #= apple(i)
        x = apple(i)
        x.printname()
        #exec("%s.init(%d)"% (x, i))
        #exec("%s.apple.printname()" %(x))
        i +=1
    """
    my_vars = {}
    for i in range(numberofclasses):
        var_name = "var%d" % i
        my_vars[var_name] = i

    print(my_vars["var2"])
        
    

    
if __name__ == '__main__':
    main()