
def maketxt():
    try:
        h, l = map(int, input("输入行数 与 列数 (例如: 3 4): ").split())
    except ValueError:
        print("格式错误")
        return None
    matrix = []
    for i in range(h):
        while True:
            line = input(f"第{i+1}行: ")
            part = line.strip().split()
            if len(part) != l:
                print(f"需要正好 {l} 个元素, 重新输入")
                continue
            try:
                row = [float(x) for x in part]
            except ValueError:
                print("包含非数字, 重输")
                continue
            matrix.append(row)
            break
    print("矩阵内容：")
    for row in matrix:
        print(" ".join(f"{x}" for x in row) )
    o = input("确认保存? (y/n): ").strip().lower()
    if o == "y":
        filename = input("输出文件名(回车默认 data.txt): ").strip() or "data"
        try:
            with open(f"{filename}.txt", "w", encoding="utf-8") as f:
                for row in matrix:
                    f.write(" ".join(f"{x}" for x in row) + "\n")
            print(f"已保存到 {filename}.txt")
        except OSError as e:
            print(f"保存失败: {e}")
    elif o == "n":
        print("已取消保存")

def readtxt(filename):
    try:
        with open(filename, "r", encoding="utf-8") as f:
            lines = f.readlines()
    except OSError as e:
        print(f"读取失败: {e}")
        return None
    matrix = []
    for line in lines:
        parts = line.strip().split()
        try:
            row = [float(x) for x in parts]
        except ValueError:
            print("文件内容包含非数字")
            return None
        matrix.append(row)
    width = len(matrix[0])
    for row in matrix:
        if len(row) != width:
            print("文件内容行长度不一致")
            return None
    return matrix

def main():
    print("请选择读取文件或是编写文件:")
    s = input("1=读取 2=编写:").strip()
    if s == "1":
        file = input("输入文件名(回车默认 data.txt): ").strip() or "data"
        readtxt(filename = "data.txt")
        matrix = readtxt(f"{file}.txt")
        if matrix is not None:
            print("矩阵内容：")
            for row in matrix:
                print(" ".join(f"{x}" for x in row) )
    elif s == "2":
        maketxt()
    else:
        print("无效选项")
        return

if __name__ == '__main__':
    main()
