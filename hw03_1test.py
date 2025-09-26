def get_matrix(path='matrix.txt'):
    # 录入方阵（n×n），并写入文件
    try:
        n = int(input("请输入方阵的阶数 n: ").strip())
    except:
        print("输入错误：需要整数。")
        return
    if n <= 0:
        print("n 必须为正整数。")
        return

    matrix = []
    for i in range(n):
        while True:
            line = input("第 {} 行（用空格分隔输入 {} 个数）: ".format(i + 1, n)).strip()
            parts = line.split()
            if len(parts) != n:
                print("需要正好 {} 个元素，请重试。".format(n))
                continue
            ok = True
            row = []
            for x in parts:
                try:
                    row.append(float(x))
                except:
                    ok = False
                    break
            if not ok:
                print("包含非数字，请重试。")
                continue
            matrix.append(row)
            break

    try:
        f = open(path, "w", encoding="utf-8")
        for row in matrix:
            f.write(" ".join(str(x) for x in row) + "\n")
        f.close()
        print("矩阵已保存到 '{}'".format(path))
    except Exception as e:
        print("写入失败: {}".format(e))

def get_vector(path='vvector.txt'):
    # 录入列向量 b，可水平(h)或垂直(l)输入，统一写到同一文件
    kind = input("向量输入方式 水平(h)/垂直(l): ").strip().lower()
    vals = []

    if kind == "h":
        line = input("请输入水平向量的全部元素（以空格分隔）: ").strip()
        parts = line.split()
        for x in parts:
            try:
                vals.append(float(x))
            except:
                print("包含非数字。")
                return
    elif kind == "l":
        print("逐行输入垂直向量元素，直接回车结束：")
        while True:
            line = input().strip()
            if line == "":
                break
            try:
                vals.append(float(line))
            except:
                print("包含非数字，请重试。")
                return
    else:
        print("无效选项。")
        return

    try:
        f = open(path, "w", encoding="utf-8")
        f.write(" ".join(str(x) for x in vals) + "\n")
        f.close()
        print("向量已保存到 '{}'".format(path))
    except Exception as e:
        print("写入失败: {}".format(e))

def read_matrix(path='matrix.txt'):
    # 从文件读取矩阵，要求为方阵
    rows = []
    try:
        f = open(path, "r", encoding="utf-8")
        for line in f:
            line = line.strip()
            if line == "":
                continue
            parts = line.split()
            row = []
            for x in parts:
                row.append(float(x))
            rows.append(row)
        f.close()
    except Exception as e:
        raise RuntimeError("读取矩阵失败: {}".format(e))

    if len(rows) == 0:
        raise ValueError("矩阵为空。")

    n = len(rows)
    for r in rows:
        if len(r) != n:
            raise ValueError("矩阵必须为方阵。")
    return rows

def read_vector(path='vvector.txt'):
    # 从文件读取向量（空白分隔的一行或多行均可）
    vals = []
    try:
        f = open(path, "r", encoding="utf-8")
        for line in f:
            line = line.strip()
            if line == "":
                continue
            parts = line.split()
            for x in parts:
                vals.append(float(x))
        f.close()
    except Exception as e:
        raise RuntimeError("读取向量失败: {}".format(e))

    if len(vals) == 0:
        raise ValueError("向量为空。")
    return vals

def gauss_solve(A, b, tol=1e-12):
    # 使用部分选主元的高斯消元法解 Ax=b
    n = len(A)
    if len(b) != n:
        raise ValueError("向量维度与矩阵不匹配。")

    # 构造增广矩阵 [A | b]
    M = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(float(A[i][j]))
        row.append(float(b[i]))
        M.append(row)

    # 前向消元（选取当前列绝对值最大的主元并交换行）
    for i in range(n):
        # 选主元行
        piv = i
        max_abs = abs(M[i][i])
        r = i + 1
        while r < n:
            if abs(M[r][i]) > max_abs:
                max_abs = abs(M[r][i])
                piv = r
            r += 1
        if max_abs < tol:
            raise ValueError("矩阵奇异或接近奇异。")

        # 交换到当前行
        if piv != i:
            tmp = M[i]
            M[i] = M[piv]
            M[piv] = tmp

        # 消元
        r = i + 1
        while r < n:
            if abs(M[i][i]) < tol:
                raise ValueError("消元阶段遇到零主元。")
            factor = M[r][i] / M[i][i]
            if abs(factor) >= tol:
                c = i
                while c <= n:
                    M[r][c] = M[r][c] - factor * M[i][c]
                    c += 1
            r += 1

    # 回代
    x = []
    for _ in range(n):
        x.append(0.0)

    i = n - 1
    while i >= 0:
        s = 0.0
        j = i + 1
        while j < n:
            s = s + M[i][j] * x[j]
            j += 1
        if abs(M[i][i]) < tol:
            raise ValueError("回代阶段遇到零主元。")
        x[i] = (M[i][n] - s) / M[i][i]
        i -= 1

    return x

def save_vector(path, vec):
    # 将解向量写入文件（空格分隔）
    try:
        f = open(path, "w", encoding="utf-8")
        f.write(" ".join(str(v) for v in vec) + "\n")
        f.close()
        print("解向量已保存到 '{}'".format(path))
    except Exception as e:
        print("写入失败: {}".format(e))

def gauss_elimination():
    # 从文件读取 -> 求解 -> 输出与保存
    try:
        A = read_matrix('matrix.txt')
        b = read_vector('vvector.txt')
        x = gauss_solve(A, b)
        print("解向量 x:")
        for i in x:
            print(f"{i:.5f}")#保留两位小数
        save_vector('result1.txt', x)
    except Exception as e:
        print("错误:", e)

def doolittle_lu(A, tol=1e-12):
    # Doolittle 分解：A = L * U，L 为单位下三角，U 为上三角
    n = len(A)
    for row in A:
        if len(row) != n:
            raise ValueError("矩阵必须为方阵。")

    L = []
    U = []
    i = 0
    while i < n:
        L.append([0.0] * n)
        U.append([0.0] * n)
        i += 1

    i = 0
    while i < n:
        L[i][i] = 1.0
        i += 1

    j = 0
    while j < n:
        # 计算 U 的第 j 列的前 j+1 行
        i = 0
        while i <= j:
            s = 0.0
            k = 0
            while k < i:
                s = s + L[i][k] * U[k][j]
                k += 1
            U[i][j] = float(A[i][j]) - s
            i += 1

        # 检查主元
        if abs(U[j][j]) < tol:
            raise ValueError("分解失败：遇到零主元（可尝试改用高斯消元或带主元的 LU）。")

        # 计算 L 的第 j 列的下三角部分
        i = j + 1
        while i < n:
            s = 0.0
            k = 0
            while k < j:
                s = s + L[i][k] * U[k][j]
                k += 1
            L[i][j] = (float(A[i][j]) - s) / U[j][j]
            i += 1

        j += 1

    return L, U

def forward_subst(L, b):
    # 解 Ly = b，其中 L 为单位下三角
    n = len(L)
    if len(b) != n:
        raise ValueError("向量维度与矩阵不匹配。")

    y = []
    i = 0
    while i < n:
        y.append(0.0)
        i += 1

    i = 0
    while i < n:
        s = 0.0
        j = 0
        while j < i:
            s = s + L[i][j] * y[j]
            j += 1
        y[i] = b[i] - s  # 因为 L 的对角为 1
        i += 1
    return y

def back_subst(U, y, tol=1e-12):
    # 解 Ux = y，其中 U 为上三角
    n = len(U)
    if len(y) != n:
        raise ValueError("向量维度与矩阵不匹配。")

    x = []
    i = 0
    while i < n:
        x.append(0.0)
        i += 1

    i = n - 1
    while i >= 0:
        s = 0.0
        j = i + 1
        while j < n:
            s = s + U[i][j] * x[j]
            j += 1
        if abs(U[i][i]) < tol:
            raise ValueError("回代阶段遇到零主元。")
        x[i] = (y[i] - s) / U[i][i]
        i -= 1
    return x

def doolittle_solve(A, b, tol=1e-12):
    # 用 Doolittle 分解解 Ax=b
    L, U = doolittle_lu(A, tol)
    y = forward_subst(L, b)
    x = back_subst(U, y, tol)
    return x

def doolittle_elimination():
    # 从 'matrix.txt' 与 'vvector.txt' 读取，求解并输出到 'result1.txt'
    try:
        A = read_matrix('matrix.txt')
        b = read_vector('vvector.txt')
        x = doolittle_solve(A, b)
        print("解向量 x:")
        for v in x:
            print(f"{v:.5f}")
        save_vector('result2.txt', x)
    except Exception as e:
        print("错误:", e)

def gauss_seidel(A, b, x0=None, tol=1e-10, max_iter=10000):
    # 用 Gauss-Seidel 迭代法解 Ax=b，返回 (x, iters)
    n = len(A)
    if len(b) != n:
        raise ValueError("向量维度与矩阵不匹配。")

    # 初始向量
    if x0 is None:
        x = []
        i = 0
        while i < n:
            x.append(0.0)
            i += 1
    else:
        if len(x0) != n:
            raise ValueError("初始向量维度与矩阵不匹配。")
        x = []
        i = 0
        while i < n:
            x.append(float(x0[i]))
            i += 1

    k = 0
    while k < max_iter:
        x_old = x[:]  # 保存上一轮用于计算未更新分量
        max_diff = 0.0

        i = 0
        while i < n:
            if abs(A[i][i]) < 1e-15:
                raise ValueError("对角元为零，无法进行 Gauss-Seidel 迭代。")
            s = 0.0
            j = 0
            while j < n:
                if j != i:
                    # j<i 用已更新的 x[j]；j>i 用上一轮的 x_old[j]
                    s = s + A[i][j] * (x[j] if j < i else x_old[j])
                j += 1
            x_new_i = (b[i] - s) / A[i][i]
            d = abs(x_new_i - x[i])
            if d > max_diff:
                max_diff = d
            x[i] = x_new_i
            i += 1

        if max_diff < tol:
            return x, (k + 1)
        k += 1

    # 未在 max_iter 内达到容差，返回当前解并提示
    # 可改为 raise ValueError("Gauss-Seidel 未收敛") 按需处理
    return x, k

def gauss_seidel_elimination():
    # 从 'matrix.txt' 与 'vvector.txt' 读取，Gauss-Seidel 求解并输出到 'result1.txt'
    try:
        A = read_matrix('matrix.txt')
        b = read_vector('vvector.txt')
        x, iters = gauss_seidel(A, b, x0=None, tol=1e-10, max_iter=10000)
        print("迭代次数:", iters)
        print("解向量 x:")
        i = 0
        while i < len(x):
            print(f"{x[i]:.5f}")
            i += 1
        save_vector('result3.txt', x)
    except Exception as e:
        print("错误:", e)

def sor(A, b, w=1.2, x0=None, tol=1e-10, max_iter=10000):
    # SOR 超松弛法解 Ax=b，返回 (x, iters)
    n = len(A)
    if len(b) != n:
        raise ValueError("向量维度与矩阵不匹配。")
    if w <= 0.0 or w >= 2.0:
        raise ValueError("松弛因子 w 需满足 0<w<2。")

    # 初始解
    if x0 is None:
        x = []
        i = 0
        while i < n:
            x.append(0.0)
            i += 1
    else:
        if len(x0) != n:
            raise ValueError("初始向量维度与矩阵不匹配。")
        x = []
        i = 0
        while i < n:
            x.append(float(x0[i]))
            i += 1

    k = 0
    while k < max_iter:
        x_old = x[:]  # 保存上一轮
        max_diff = 0.0

        i = 0
        while i < n:
            aii = A[i][i]
            if abs(aii) < 1e-15:
                raise ValueError("对角元为零，无法进行 SOR 迭代。")

            # s1 = Σ_{j<i} a_ij * x_j(新)
            s1 = 0.0
            j = 0
            while j < i:
                s1 = s1 + A[i][j] * x[j]
                j += 1

            # s2 = Σ_{j>i} a_ij * x_j(旧)
            s2 = 0.0
            j = i + 1
            while j < n:
                s2 = s2 + A[i][j] * x_old[j]
                j += 1

            x_gs = (b[i] - s1 - s2) / aii
            x_new_i = (1.0 - w) * x[i] + w * x_gs

            d = abs(x_new_i - x[i])
            if d > max_diff:
                max_diff = d
            x[i] = x_new_i
            i += 1

        if max_diff < tol:
            return x, (k + 1)
        k += 1

    return x, k  # 未在 max_iter 内达到容差

def sor_elimination():
    # 从 'matrix.txt' 与 'vvector.txt' 读取，SOR 求解并输出到 'result4.txt'
    try:
        A = read_matrix('matrix.txt')
        b = read_vector('vvector.txt')
        # 可调整 w 来改善收敛速度（常在 1.0~1.9 之间试探）
        w = 1.2
        x, iters = sor(A, b, w=w, x0=None, tol=1e-10, max_iter=10000)
        print("SOR 松弛因子 w:", w)
        print("迭代次数:", iters)
        print("解向量 x:")
        i = 0
        while i < len(x):
            print(f"{x[i]:.5f}")
            i += 1
        save_vector('result4.txt', x)
    except Exception as e:
        print("错误:", e)

if __name__ == "__main__":
    # 如需交互创建数据可先运行：
    # get_matrix()
    # get_vector()
    gauss_elimination()
    doolittle_elimination()
    gauss_seidel_elimination()
    sor_elimination()