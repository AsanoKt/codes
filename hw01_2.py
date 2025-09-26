def getmatrix():
    print("输入矩阵：")
    h, l = map(int,input("行数和列数:").split())
    matrix = []
    for i in range(h):
        while True:
            line = input(f"第{i+1}行:").strip().split()
            if len(line) != l:
                print(f"需要正好 {l} 个元素, 重新输入")
                continue
            try:
                row = [float(x) for x in line]
            except ValueError:
                print("包含非数字, 重输")
                continue
            matrix.append(row)
            break
    print("矩阵内容：")
    for row in matrix:
        print(" ".join(f"{x}" for x in row))
    '''print(matrix)'''
    return matrix

def getvector():
    print("输入向量：")
    vtype = input("行向量/(默认)列向量 (h/l):").strip().lower() or "l"
    vector = []
    if vtype == "h":
        while True:
            line = input("用空格分隔的数字:").strip().split()
            try:
                vector = [float(x) for x in line]
                break
            except ValueError:
                print("包含非数字, 重输")
                continue
        print("向量内容：")
        print(" ".join(f"{x}" for x in vector))
    elif vtype == "l":
        print("逐行输入数字, 空行结束:")
        while True:
            line = input().strip()
            if line == "":
                break
            try:
                num = float(line)
                vector.append(num)
            except ValueError:
                print("包含非数字, 重输")
                continue
        print("向量内容：")
        for x in vector:
            print(f"{x}")
    '''print(vector)'''
    return vector

def mat_vec(m, v):
    if len(m[0]) != len(v):
        print("矩阵列数与向量长度不匹配")
        return None
    result = []
    for row in m:
        s = sum(row[i] * v[i] for i in range(len(v)))
        result.append(s)
    return result

def main():
    for _ in range(2):
        if _ == 0:
            m = getmatrix()
            continue
        elif _ == 1:
            v = getvector()
            continue
    mat_vec(m, v)
    r = mat_vec(m, v)
    if r is not None:
        print("结果向量：")
        for x in r:
            print(f"{x}")

if __name__ == '__main__':
    main()