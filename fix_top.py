#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fix_topol.py

用途：修正当前目录下 topol.top 文件的格式问题：
1. 将所有 #include 语句拆分为独立行；
2. 补全 [ system ] 块，体系名称取自父目录名（即 MOF 名称）。

使用方法：
    在包含 topol.top 的目录（例如 ACN/ 或 PC/）中直接运行：
        python fix_topol.py
"""

import os
import sys
import re

def extract_include_lines(lines):
    """
    从所有行中提取每个 #include 语句，支持一行多个 #include。
    返回一个列表，每个元素为完整的 #include 行（以 #include 开头）。
    """
    includes = []
    for line in lines:
        if '#include' not in line:
            continue
        # 按 '#include' 分割，但保留 '#' 位置，我们使用正则或简单 split
        # 方法：找到所有 '#include' 出现的位置，截取子串
        # 简单方式：用 '#include' 分割，然后组装
        parts = line.split('#include')
        # parts[0] 可能为空白或前导空格，忽略
        for i in range(1, len(parts)):
            # 提取从 '#include' 开始到该行结束或下一个 '#include' 之前
            # 但我们已经分割了，只需处理 parts[i] 并加上 '#include'
            inc_part = parts[i].strip()
            if inc_part:
                # 需要提取到行尾，但可能后面有注释或空格
                # 去掉行内可能残留的其他 include 标记，但已分割，所以直接组合
                inc_line = '#include ' + inc_part.split('#include')[0].strip()
                # 确保是完整的语句（可能有引号）
                includes.append(inc_line)
    return includes

def find_molecules_index(lines):
    """
    查找 [ molecules ] 所在的行索引，若不存在返回 -1。
    """
    for i, line in enumerate(lines):
        if line.strip().lower().startswith('[ molecules ]'):
            return i
    return -1

def main():
    # 当前目录
    cwd = os.getcwd()
    topol_path = os.path.join(cwd, 'topol.top')
    if not os.path.isfile(topol_path):
        print(f"错误：当前目录下找不到 topol.top 文件：{cwd}")
        sys.exit(1)

    # 读取原文件
    with open(topol_path, 'r') as f:
        orig_lines = f.readlines()

    # 提取所有 include 语句
    include_lines = extract_include_lines(orig_lines)
    if not include_lines:
        print("警告：未检测到任何 #include 语句，可能文件格式异常。")

    # 定位 [ molecules ] 位置
    mol_idx = find_molecules_index(orig_lines)
    if mol_idx == -1:
        print("错误：未能找到 [ molecules ] 部分。")
        sys.exit(1)

    # 保存 [ molecules ] 及其之后的所有行
    molecules_rest = orig_lines[mol_idx:]

    # 体系名称：取父目录名作为 MOF 名称
    parent_dir = os.path.basename(os.path.dirname(cwd))
    # 如果父目录名是 '.' 或空，则用当前目录名
    if parent_dir in ('', '.'):
        parent_dir = os.path.basename(cwd)
    system_name = parent_dir

    # 构建新内容
    new_lines = []
    # 写入所有 include
    for inc in include_lines:
        new_lines.append(inc + '\n')
    new_lines.append('\n')
    # [ system ] 块
    new_lines.append('[ system ]\n')
    new_lines.append(system_name + '\n')
    new_lines.append('\n')
    # 写入 [ molecules ] 及其后内容
    new_lines.extend(molecules_rest)

    # 写入临时文件
    tmp_path = os.path.join(cwd, 'topol_new')
    with open(tmp_path, 'w') as f:
        f.writelines(new_lines)

    # 替换原文件
    os.remove(topol_path)
    os.rename(tmp_path, topol_path)

    print(f"修正完成：{topol_path} 已更新")
    print(f"体系名称设置为：{system_name}")
    print(f"共处理了 {len(include_lines)} 个 #include 语句。")

if __name__ == '__main__':
    main()
