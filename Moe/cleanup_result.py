

def main():
    filename = "initial_test_results_12_09/times_papso.dat"

    with open(filename, "r") as f:
        lines = f.readlines()
    with open(filename, "w") as f:
        for line in lines:
            if len(line.strip("\n").split(" ")) == 4:
                f.write(line) 

if __name__ == '__main__':
    main()
