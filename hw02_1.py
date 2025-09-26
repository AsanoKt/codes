import time

def J_compute():#雅各比迭代法
    xk = float(input("initial x="))#获取输入内容
    start_time = time.time()
    xk = abs(xk)#取绝对值
    xk1 = (xk * 3.0) ** (1/3) #迭代公式
    '''print(xk1)'''
    i = 0
    while abs(xk1 - xk) >= 1e-6: #精度要求
        xk = xk1
        '''print(f"xk{xk}")'''
        xk1 = (xk * 3.0) ** (1 / 3)
        '''print(f"xk1{xk1}")'''
        i += 1#统计迭代次数
    '''xk2 = - (xk * 3.0) ** (1 / 3) #考虑负值
    while abs(xk2 - xk) >= 1e-6:
        xk = xk2
        xk2 = -(-xk * 3.0) ** (1 / 3)
        i += 1#统计迭代次数
    x0 = 0.0 #考虑0值
    xk3 = (x0 * 3.0) ** (1 / 3)
    while abs(x0 - xk3) >= 1e-6:
        x0 = xk3
        xk3 = (x0 * 3.0) ** (1 / 3)
        #print(f"xk1{xk1}")
        i += 1#统计迭代次数'''
    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.10f} seconds")
    print(f"Iteration times: {i}")#迭代次数
    return xk1 #xk3, xk2 #返回三个根

def N_compute():#牛顿法
    xk = float(input("initial x="))
    start_time = time.time()
    i = 0
    if xk >= 0:
        xk1 = xk - (xk**3 - 3*xk) / (3*xk**2 - 3)
        '''print(xk1)'''
        while abs(xk1 - xk) >= 1e-5:
            i += 1
            xk = xk1
            '''print(f"xk{xk}")'''
            xk1 = xk - (xk**3 - 3*xk) / (3*xk**2 - 3)
            '''print(f"xk1{xk1}")'''
    elif xk < 0:
        xk1 = xk - (xk**3 - 3*xk) / (3*xk**2 - 3)
        while abs(xk1 - xk) >= 1e-5:
            i += 1
            xk = xk1
            '''print(f"xk{xk}")'''
            xk1 = xk - (xk**3 - 3*xk) / (3*xk**2 - 3)
            '''print(f"xk1{xk1}")'''
    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.10f} seconds")
    print(f"Iteration times: {i}")#迭代次数
    return xk1

def post_acc(x0):#post acceleration法
    start_time = time.time()
    gx0 = x0**2
    xk = x0
    xk1 = (xk ** 3 / 3 - xk * gx0) / (1 - gx0)
    i = 0
    while abs(xk - xk1) >= 1e-5:
        xk = xk1
        xk1 = (xk ** 3 / 3 - xk * gx0) / (1 - gx0)
        i += 1#统计迭代次数

    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.10f} seconds")
    print(f"Iteration times: {i}")  # 迭代次数
    return xk1

def aitken_acc(x0):#Aitken加速法
    start_time = time.time()
    xk = x0
    xk1 = xk ** 3 / 3
    xk2 = xk1 ** 3 / 3
    i = 0
    while abs(xk2 - xk1) >= 1e-5:
        xk = xk1
        xk1 = xk ** 3 / 3
        xk2 = xk1 ** 3 / 3
        if (xk2 - 2*xk1 + xk) != 0:
            xk1 = xk - ((xk1 - xk)**2) / (xk2 - 2*xk1 + xk)
        else:
            print("除零错误，停止迭代")
            break
        i += 1

    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.10f} seconds")
    print(f"Iteration times: {i}")  # 迭代次数
    return xk1

if __name__ == "__main__":

    print("f(x) = x^3/3 - x = 0")
    m = input("method (J/N/P/A): ")
    if m == "N":
        result = N_compute()
        print(f"final x={result:.5f}")
    elif m == "J":
        r0 = J_compute()
        print(f"final x={r0:.5f}")
    elif m == "P":
        x0 = float(input("initial x="))
        result = post_acc(x0)
        print(f"final x={result:.5f}")
    elif m == "A":
        x0 = float(input("initial x="))
        result = aitken_acc(x0)
        print(f"final x={result:.5f}")
