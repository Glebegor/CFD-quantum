import qiskit as qk



class QLBM:
    def __init__(self, params):
        '''
        size: 512x512
        qreg_x: 9
        qreg_y: 9
        qreg_d: 4
        qreg: 22
        creg: 22
        
        '''
        self.qreg_x = qk.QuantumRegister(params.REG_X_SIZE)
        self.qreg_y = qk.QuantumRegister(params.REG_Y_SIZE)
        self.qreg_d = qk.QuantumRegister(params.REG_D_SIZE)
        self.creg = qk.ClassicalRegister(params.REG_X_SIZE+ params.REG_Y_SIZE + params.REG_D_SIZE)
        