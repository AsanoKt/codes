def Pmn(n,m):
    a = 1
    for i in range(n-m+1,n+1):
        a *= i
    return a

def Cmn(n,m):
    return Pmn(n,m)//Pmn(m,m)

def main():
    n = int(input("输入n:"))
    m = int(input("输入m:"))
    print(f"P({n},{m})={Pmn(n,m)}")
    print(f"C({n},{m})={Cmn(n,m)}")

if __name__ == '__main__':
    main()