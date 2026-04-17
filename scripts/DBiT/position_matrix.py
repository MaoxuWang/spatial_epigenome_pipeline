import sys
import os

def main():
    # 模拟 Perl 的 $ARGV[0]
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_file>", file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]

    # 参数设置
    pwidth = 1200
    pheight = 1200
    flank_x = 100
    flank_y = 100

    # 数据存储结构 (模拟 Perl 的 %hash)
    # Python 字典: key -> [x, y, genes_count, total_counts]
    data_hash = {}
    max_val = 0

    try:
        with open(input_file, 'r') as f:
            # <IN>; 模拟 Perl 跳过第一行 (Header)
            next(f)

            for line in f:
                # chomp; (去除末尾换行符)
                line = line.rstrip('\n')
                if not line:
                    continue

                array = line.split('\t')
                
                # key 是 array[0] (例如 "10x20")
                key = array[0]
                
                try:
                    pos = key.split('x')
                    x_pos = int(pos[0])
                    y_pos = int(pos[1])
                except ValueError:
                    # 如果坐标格式不对，跳过该行 (Perl可能会报错或产生奇怪结果，Python选择安全跳过)
                    continue

                # 初始化 hash 结构
                # 对应 Perl: $hash{$array[0]}[0], [1], [2], [3]
                if key not in data_hash:
                    # [x, y, genes, counts]
                    data_hash[key] = [x_pos, y_pos, 0, 0.0]

                for i in range(1, len(array)):
                    val_str = array[i]
                    
                    # 模拟 Perl 的数值转换行为
                    # 在 Perl 中，非数字字符串会被视为 0，这里做类似处理
                    try:
                        val = float(val_str)
                    except ValueError:
                        val = 0

                    if val > 0:
                        data_hash[key][2] += 1
                        data_hash[key][3] += val

                if data_hash[key][2] > max_val:
                    max_val = data_hash[key][2]

    except FileNotFoundError:
        print(f"Error: File {input_file} not found.", file=sys.stderr)
        sys.exit(1)

    # 开始输出 SVG XML (print $svg->xmlify())
    # 模拟 Perl SVG 模块的标准输出头
    print('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>')
    print('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">')
    print(f'<svg height="{pheight}" width="{pwidth}" xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">')

    # 背景黑色矩形
    print(f'<rect x="100" y="100" width="990" height="990" fill="black" />')

    for key in data_hash:
        item = data_hash[key]
        x = item[0]
        y = item[1]
        gene_count = item[2]

        degree = 0
        if max_val > 0:
            degree = gene_count / max_val
        
        # 坐标变换逻辑
        rect_x = flank_x + (x - 1) * 20
        rect_y = flank_y + (y - 1) * 20

        # 保留4位小数以保持XML整洁
        print(f'<rect x="{rect_x}" y="{rect_y}" width="10" height="10" opacity="{degree:.6f}" fill="red" />')

    print('</svg>')

if __name__ == "__main__":
    main()