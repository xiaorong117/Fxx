#!/usr/bin/env python3
"""Append the current implementation details to the supplied control-equation DOCX.

The script deliberately uses only Python's standard OOXML/ZIP support so the
update is reproducible in the offline production container.  It preserves every
original package part and refuses to overwrite the input or an existing output.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(text: str, *, bold: bool = False, italic: bool = False) -> str:
    props = []
    if bold:
        props.append("<w:b/>")
    if italic:
        props.append("<w:i/>")
    rpr = f"<w:rPr>{''.join(props)}</w:rPr>" if props else ""
    space = ' xml:space="preserve"' if text[:1].isspace() or text[-1:].isspace() else ""
    return f"<w:r>{rpr}<w:t{space}>{escape(text)}</w:t></w:r>"


def paragraph(text: str, style: str = "FirstParagraph",
              *, bold_prefix: str = "") -> str:
    content = ""
    if bold_prefix and text.startswith(bold_prefix):
        content = run(bold_prefix, bold=True) + run(text[len(bold_prefix):])
    else:
        content = run(text)
    return (f'<w:p><w:pPr><w:pStyle w:val="{style}"/></w:pPr>'
            f"{content}</w:p>")


def heading(text: str) -> str:
    return paragraph(text, "Heading1")


def bullet(text: str) -> str:
    return paragraph("• " + text, "Compact")


def equation(text: str) -> str:
    # Word-native OMML container.  Linear mathematical notation is intentionally
    # kept explicit and searchable rather than inserted as a raster image.
    return ("<w:p><w:pPr><w:pStyle w:val=\"BodyText\"/>"
            "<w:jc w:val=\"center\"/></w:pPr><m:oMathPara><m:oMath>"
            f"<m:r><m:rPr><m:sty m:val=\"i\"/></m:rPr><m:t>{escape(text)}</m:t></m:r>"
            "</m:oMath></m:oMathPara></w:p>")


def page_break() -> str:
    return '<w:p><w:r><w:br w:type="page"/></w:r></w:p>'


def appendix_xml() -> str:
    parts = [
        page_break(),
        heading("20. 当前实现对齐修订说明"),
        paragraph(
            "本节及第21—31节依据当前生产代码补充原文第19节之后的数值实现细节。"
            "若与原文中“建议”“可采用”等表述不一致，以本修订部分描述的实际实现为准；"
            "式（29）—（35）的全耦合控制方程本身保持不变。修订没有引入第二套自动微分"
            "或第二套物理残量。"),
        paragraph(
            "实现范围　当前程序使用同一套FADBAD++自动微分残量、稀疏Jacobian和"
            "monolithic Newton流程。活动含气控制体保留Pl、C、Pg、Sw；无气节点仅保留"
            "Pl、C。自由气体没有网络通量，气泡仍为静止、不可迁移和不可重新激活。",
            "BlockText"),
        heading("21. 相对压力、绝对压力与传质倍率"),
        paragraph(
            "求解未知量和VTK中的液相、气相压力采用相对压力。状态方程和Henry平衡必须"
            "恢复同一个背景绝对压力Pbg。毛细压差中背景压力相消。"),
        equation("P_l^abs = P_bg + P_l,    P_g^abs = P_bg + P_g.    (43)"),
        equation("n_g,i = P_g,i^abs V_i(1-S_w,i)/(ZRT).    (44)"),
        paragraph(
            "界面传质允许显式倍率Mtr，但该倍率必须同时作用于溶解组分残量和自由气体残量"
            "中的同一个通量。科研基准采用Mtr=1；Mtr≠1只能标记为加速预览。"),
        equation("J_gl,i = M_tr k_L A_gl,i [H_cp P_g,i^abs - C_i].    (45)"),
        heading("22. 总入口流量方程、浓度边界与局部回流"),
        paragraph(
            "流量入口不是逐喉指定流量。程序增加一个共享入口压力Pin作为未知量，并增加一条"
            "全局方程；各入口喉的流量由局部导流系数和压力解自然分配。出口液相压力固定。"),
        equation("R_Q = Σ_(e∈Γ_in) G_e(P_in-P_e) - Q_in,target = 0.    (46)"),
        paragraph(
            "入口正向流入使用Cin；入口局部反流使用内部上游浓度。出口正常流出使用内部浓度，"
            "出口回流使用Cbackflow；出口扩散采用零梯度。总入口流量为正并不排除个别入口喉"
            "反流，也不排除个别出口喉回流。程序分别记录入口反流和出口回流的边数及体积流率。"),
        paragraph(
            "气泡缩小时液体体积增加，因此瞬时入口和出口流量一般不相等。其差额必须与内部"
            "液体体积变化及active-set体积账本共同满足体积守恒，不能强制Qin=Qout。"),
        heading("23. 零距离cluster、massless junction与压力规范"),
        paragraph(
            "零半长度连接不被替换成人为正长度。程序用union-find构造零距离cluster；cluster"
            "共享一个Pl和C，并组装一个聚合Rl、Rd。每个正体积含气成员仍保留独立Pg、Sw、"
            "EOS库存和界面面积。"),
        paragraph(
            "零体积且Sw=1的控制体作为massless junction存在，不含积累项。零体积却含气、"
            "负体积、负长度仍为硬错误。与压力边界隔离且完全充液的分量具有压力常数零空间，"
            "程序用上一接受状态的一个液相压力规范方程替换冗余Rl；该规范不改变流量和守恒。"),
        heading("24. 自动微分Newton流程与双收敛门槛"),
        paragraph(
            "每个Newton迭代都重新评价上风方向、残量和Jacobian数值；只允许线性后端复用"
            "固定稀疏结构，不冻结Jacobian。fraction-to-boundary保持Pg绝对压力为正、C非负且"
            "Sw位于允许范围，Armijo回溯要求缩放残量下降。"),
        equation("||R_scaled||_∞ ≤ ε_R.    (47)"),
        equation("||ΔU_scaled||_∞ ≤ ε_U [1+||U_scaled||_∞].    (48)"),
        paragraph(
            "只有式（47）与式（48）同时满足，且状态边界、非负自由气体、摩尔守恒和液体"
            "体积守恒全部通过，时间步才接受。线性API返回成功不能替代CPU端真实||Ax-b||检查。"
            "Eigen、AMGX和PARDISO只替换线性求解环节，不改变控制方程。请求不可用后端时明确"
            "失败，不静默回退。"),
        heading("25. 自适应时间步控制"),
        paragraph(
            "自适应控制器扩展原有失败重试路径，不与其并列建立第二套控制器。设本步接受步长"
            "为Δt、Newton迭代数为k，默认控制为："),
        bullet("k≤3：候选下一步为min(1.25Δt, Δtmax)，连续增长不超过3次；"),
        bullet("4≤k≤7：保持Δt；"),
        bullet("k≥8：候选下一步为max(0.8Δt, Δtmin)；"),
        bullet("Newton、线性、状态或守恒失败：拒绝本步，以0.5Δt重算；"),
        bullet("达到Δtmin仍不能接受：明确失败，禁止静默接受。"),
        paragraph(
            "气体时间尺度不再使用全网最小ng/|dng/dt|。程序排除已进入消失邻域及自由气体摩尔"
            "占比低于配置门槛的泡体，对其余有效泡体的τg=ng/|dng/dt|取配置分位数（当前默认"
            "中位数），再施加gas_timescale_fraction。这样单个微小气泡不会永久控制全局步长。"),
        paragraph(
            "近零局部浓度不允许通过C/|dC/dt|压缩全局步长。浓度变化限制默认关闭；若开启，"
            "它使用明确物理尺度Cref=Hcp Pbg约束|ΔCi|≤fc Cref，并输出限制节点和候选步长。"),
        heading("26. 单气泡精确消失事件定位"),
        paragraph(
            "完成隐式候选步后，程序检查每个原活动泡体的新自由气体摩尔数和Sg。预测负摩尔数"
            "会直接拒绝，不能靠截断掩盖。若跨越ng,off或Sg,off且不符合批量退休条件，则根据"
            "旧、新状态线性估计最早阈值穿越比例fi："),
        equation("f_i = min{(n_g^n-n_off)/(n_g^n-n_g*), (S_g^n-S_g,off)/(S_g^n-S_g*)}.    (49)"),
        paragraph(
            "程序用略跨过阈值的事件步Δtevent=clamp(1.0001 fi Δt, Δtmin, Δt)重新求解。只有事件"
            "步通过Newton、状态和守恒后，才将该泡体切换为永久无气状态。事件后下一请求步长"
            "从事件前请求步长的安全比例恢复，而不是从极小事件步长进行数千次几何增长。"),
        heading("27. 守恒型微小气泡批量退休"),
        paragraph(
            "为避免大量数值上可忽略的气泡逐一定位并卡住全局时间，程序提供可选批量退休。"
            "设候选泡体旧摩尔数为ng,i，当前全局自由气体为Ng。只有同时满足下列硬门槛才批量："),
        equation("n_g,i/N_g ≤ f_ind,    Σ_(i∈B)n_g,i/N_g ≤ f_batch.    (50)"),
        paragraph(
            "当前已验证配置采用find=10^-10、fbatch=10^-8。任一单泡超限、候选总量超限或出现"
            "负摩尔预测时不走批量路径，而保留精确定位。阈值可配置，但不能为了速度任意放宽。"),
        paragraph(
            "批量退休不是删除气体：接受步末仍为正的自由气体逐泡加入所属cluster的溶解库存，"
            "相应气体体积进入active-set液体体积账本，然后重新执行原摩尔和体积守恒硬门槛。"),
        equation("N_d,cluster ← N_d,cluster + Σ_(i∈B)n_g,i*,    V_ledger ← V_ledger + Σ_(i∈B)V_i S_g,i*.    (51)"),
        paragraph(
            "事件表中的estimated_crossing_time_s由接受步内插值得到；它适合记录时空分布，但"
            "不是在该时刻重新求解完整耦合系统得到的严格事件时间。论文若分析单泡精确消失"
            "时间，应另设关闭批量路径的代表性对照。"),
        heading("28. active-set切换与守恒审计"),
        paragraph(
            "事件接受后设置active_gas=false、permanently_inactive=true、Sw=1、Pg=NaN（表示该"
            "自由度不存在），并在下一步重新建立活动自由度和稀疏图。无气节点继续求解Pl、C。"),
        equation("N_0 + N_in - N_out + N_ledger = N_free + N_dissolved + E_N.    (52)"),
        equation("V_l,0 + V_l,in - V_l,out + V_ledger = V_l(t) + E_V.    (53)"),
        paragraph(
            "入口和出口累计量均以正的gross量保存；净补液另行记录。epsilon_N_adjusted和"
            "epsilon_V_adjusted使用最终收敛状态、实际边界通量及active-set账本计算。方程残量"
            "通过、线性收敛或程序正常退出均不能替代独立守恒审计。"),
        heading("29. 科学统计、区段与终止条件"),
        paragraph(
            "总体与主贯通域饱和度只统计真实、正体积控制体；排除虚拟储层和massless junction，"
            "每个生产控制体体积只计一次。主贯通域定义为入口—出口连通的主要分量。"),
        equation("S̄_g,main = [Σ_(i∈Ω_main,real)V_i(1-S_w,i)]/[Σ_(i∈Ω_main,real)V_i].    (54)"),
        paragraph(
            "支持S̄g,main≤Sg,target和Nfree(t)≤fres Nfree(0)两类科学终止，并保留最大物理时间、"
            "最大墙钟、磁盘保护、安全STOP及数值失败原因。最大物理时间先达到属于右删失，"
            "不能报告为气泡已完全溶解。"),
        paragraph(
            "默认沿归一化流向坐标划分20段。孔按孔心唯一归段，喉按中点归段；每段记录真实"
            "孔隙体积、Sg、气体体积、自由/溶解摩尔数、界面面积、累计传质、活动气泡数、"
            "消失事件及跨段通量。所有分段量之和必须与对应全局量一致。"),
        heading("30. checkpoint、输出与可复现性"),
        paragraph(
            "checkpoint保存状态、时间、步号、next_dt、连续增长计数、入口共享压力、累计边界"
            "通量、相间传质、active-set账本、逐节点累计传质、事件计数、批量退休配置和不可逆"
            "状态。restart校验节点数、背景压力、传质倍率、延续配置签名及PNM哈希；不兼容时"
            "明确失败。restart后TSV不得重复时间、丢失累计量或重置分段积分。"),
        paragraph(
            "diagnostics.tsv逐接受步记录Newton、线性残量、四残量、dt原因、库存、流量、事件"
            "和守恒；global_history.tsv用于全局曲线；segment_history.tsv用于时空统计；"
            "gas_disappearance_events.tsv逐泡记录；balance_history.tsv保存独立守恒。VTP/PVD"
            "输出相对压力、Sw/Sg、C、active状态、气体/溶解摩尔数、边界标记和传质倍率。"),
        heading("31. 验证边界与论文报告要求"),
        paragraph(
            "数值验收应依次覆盖小线性系统、小网络、固定/自适应dt对照、精确/批量事件对照、"
            "连续/restart一致性和代表性全网络。每个后端必须满足相同Newton双门槛、CPU真实"
            "线性残量及守恒门槛。性能必须报告端到端时间，不能只报告solve_seconds。"),
        paragraph("论文至少应区分下列三类内容：", "FirstParagraph"),
        bullet("控制方程和严格守恒耦合：当前代码已实现并有回归/机器输出支持；"),
        bullet("数值策略：自适应dt、事件定位、批量退休和后端属于离散求解方法；"),
        bullet("物理闭合：Hcp、kL、界面面积、Qin半径/b、接触角、相对导流、扩散面积/曲折度、"
               "静止不可再活化气泡及选定Pe/Ca仍需文献或实验标定。"),
        paragraph(
            "因此，守恒型批量退休可以用于大规模级配比较中的数值加速，但其阈值、累计退休"
            "质量比例、最大epsilon_N/epsilon_V以及代表性精确事件对照必须在论文中披露。",
            "BlockText"),
    ]
    return "".join(parts)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    source = args.input.resolve()
    output = args.output.resolve()
    audit = args.audit.resolve()
    if source == output or output.exists():
        raise SystemExit("refusing to overwrite input or existing output")
    if audit.exists():
        raise SystemExit("refusing to overwrite existing audit")
    with ZipFile(source) as zin:
        names = zin.namelist()
        document = zin.read("word/document.xml").decode("utf-8")
        if "20. 当前实现对齐修订说明" in document:
            raise RuntimeError("input already contains the implementation appendix")
        document = document.replace(
            "全耦合控制方程、离散形式与数值实现",
            "全耦合控制方程、离散形式与当前数值实现（实现对齐修订版）",
            1)
        marker = "<w:sectPr"
        position = document.rfind(marker)
        if position < 0:
            raise RuntimeError("DOCX document.xml has no final sectPr")
        document = document[:position] + appendix_xml() + document[position:]
        core = zin.read("docProps/core.xml").decode("utf-8")
        core = re.sub(
            r"(<dc:title>).*?(</dc:title>)",
            r"\1静止气泡溶解—全耦合控制方程与当前数值实现（修订版）\2",
            core, count=1)
        now = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace(
            "+00:00", "Z")
        core = re.sub(
            r"(<dcterms:modified[^>]*>).*?(</dcterms:modified>)",
            rf"\g<1>{now}\g<2>", core, count=1)
        output.parent.mkdir(parents=True, exist_ok=True)
        with ZipFile(output, "w", ZIP_DEFLATED) as zout:
            for name in names:
                data = zin.read(name)
                if name == "word/document.xml":
                    data = document.encode("utf-8")
                elif name == "docProps/core.xml":
                    data = core.encode("utf-8")
                zout.writestr(name, data)
    # Structural readback independent of the generation strings.
    with ZipFile(output) as check:
        required_parts = {"[Content_Types].xml", "word/document.xml",
                          "word/styles.xml", "word/settings.xml",
                          "_rels/.rels", "docProps/core.xml"}
        missing = sorted(required_parts - set(check.namelist()))
        final_xml_bytes = check.read("word/document.xml")
        ET.fromstring(final_xml_bytes)
        ET.fromstring(check.read("docProps/core.xml"))
        ET.fromstring(check.read("[Content_Types].xml"))
        final_xml = final_xml_bytes.decode("utf-8")
        if missing or final_xml.count("<w:sectPr") != 1:
            raise RuntimeError(f"invalid output DOCX package; missing={missing}")
        for required in ("20. 当前实现对齐修订说明", "27. 守恒型微小气泡批量退休",
                         "estimated_crossing_time_s", "epsilon_N_adjusted",
                         "31. 验证边界与论文报告要求"):
            if required not in final_xml:
                raise RuntimeError(f"DOCX readback missing section: {required}")
    result = {
        "format": "bubble_control_equations_docx_update_audit_v1",
        "status": "passed",
        "input": str(source),
        "input_sha256": sha256(source),
        "output": str(output),
        "output_sha256": sha256(output),
        "method": "preserve original DOCX package and append native OOXML paragraphs/math",
        "original_package_parts": len(names),
        "new_sections": list(range(20, 32)),
        "equation_numbers_added": list(range(43, 55)),
        "original_overwritten": False,
        "updated_utc": now,
    }
    audit.parent.mkdir(parents=True, exist_ok=True)
    audit.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
