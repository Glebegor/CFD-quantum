import params
import qlbm



class APP:
    def __init__(self):
        self.params = params
        self.qlbm = qlbm.QLBM(self.params)
                
    def printIntro(self):
        print(f"##" + ("-"*40) + "##")
        print(f"Creators: Hlib Arseniuk and Jiri Posavad")
        print(f"Project: QLBM for CFD")
        print(f"Size of cells X*X area: " + str(self.params.CELLS_SIZE))
        print(f"Register: " + str(self.params.REG_X_SIZE + self.params.REG_Y_SIZE + self.params.REG_D_SIZE))
        print(f"Directions: ")
        print(f"6  2  5")
        print(f" \ | / ")
        print(f"3  0  1")
        print(f" / | \ ")
        print(f"7  4  8")
        print(f"##" + ("-"*40) + "##")
    
    def Run(self):
        self.printIntro()
    
    
    
if __name__ == "__main__":
    app = APP()
    app.Run()