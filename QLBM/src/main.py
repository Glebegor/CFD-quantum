import argparse

import params
import params_test
import qlbm
import qiskit_ibm_runtime as qk_runtime
import credentials



class APP:
    def __init__(self, selected_params=params, time_series=False):
        self.params = selected_params
        self.qlbm = qlbm.QLBM(self.params)
        self.service = None
        self.time_series = time_series
                
    def printIntro(self):
        print(f"##" + ("-"*40) + "##")
        print(f"Creators: Hlib Arseniuk and Jiri Posavad")
        print(f"Project: QLBM for CFD")
        print(f"Mode: " + ("LOCAL" if self.params is params_test else "FULL"))
        print(f"Size of cells X*X area: " + str(self.params.CELLS_SIZE))
        print(f"Register: " + str(self.params.REG_X_SIZE + self.params.REG_Y_SIZE + self.params.REG_D_SIZE))
        print(f"Directions: ")
        print(f"6  2  5")
        print(" \\ | / ")
        print(f"3  0  1")
        print(" / | \\ ")
        print(f"7  4  8")
        print(f"Type of parameters: {self.params.TYPE}")
        print(
            "Time series: "
            f"{getattr(self.params, 'SNAPSHOT_COUNT', 40)} snapshots × "
            f"{getattr(self.params, 'TIME_STEP', 0.1)} s; "
            f"{getattr(self.params, 'SHOTS', 1000)} measurement shots/snapshot"
        )
        print(f"##" + ("-"*40) + "##")
    
    def setupIBMApi(self):
        self.service = qk_runtime.QiskitRuntimeService(
            token=credentials.TOKEN,
            instance=credentials.CRN,
            channel="ibm_quantum_platform"
        )
    
    def Run(self):
        self.printIntro()
        if self.params.TYPE == "prod":
            self.setupIBMApi()
        if self.time_series:
            self.qlbm.RunTimeSeries(self.service)
        else:
            self.qlbm.Run(self.service)
    
    
    
def ParseArguments():
    parser = argparse.ArgumentParser(description="Run the QLBM application")
    parser.add_argument(
        "-local",
        "--local",
        action="store_true",
        help="use the small local 4x4 configuration",
    )
    parser.add_argument(
        "--time-series",
        action="store_true",
        help="run 40 snapshots at 0.1-second intervals",
    )
    return parser.parse_args()


if __name__ == "__main__":
    
    
    arguments = ParseArguments()
    app = APP(
        params_test if arguments.local else params,
        time_series=arguments.time_series,
    )
    app.Run()
