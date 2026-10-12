# From the Sound Barrier to Relativistic Inertia

**Retarded Fields and an Effective-Continuum Conjecture**

**An exploratory physics note by Mingxiang Shen (沈茗祥), Independent Researcher.**

[Read the paper (PDF, 17 pages)](output/pdf/wave_speed_inertia.pdf) · [LaTeX source](manuscript.tex) · [Verification script](verify_calculations.py) · [Recorded checks](verification_results.txt)

## 从音障到光速极限：一个业余物理联想

这个想法来自一个简单的问题：飞机接近声速时，有限的声波传播速度会改变流场并产生跨声速效应；物体接近光速时，是否也可能与有限速度传播的场的动态响应存在某种联系？普朗特—格劳厄特修正因子与洛伦兹因子相似的数学形式，是这一联想的起点。

这份笔记沿着这个问题做了一些明确的模型计算，并进一步猜想：时空几何、粒子和惯性是否可能共同来自某种底层连续系统的有效描述。

作者是一名计算机专业学生，这项探索来自业余思考，并在与 OpenAI GPT 的讨论中形成。GPT 实质参与了推导、文献查找、组织和写作。这里公开保存思考过程、稿件和可复现检查，欢迎有兴趣的读者指出问题、补充相关研究。

### 如何理解这份笔记

- **模型中的推导：** 普朗特—格劳厄特修正、相对论纵向和横向响应、匀速源的延迟势、线性化引力、束缚能、四极辐射、声学有效度规，以及 sine-Gordon 孤子的能量和动量。
- **目前的结论：** 文中的标准弱场计算不支持把“前方引力波堆积”直接作为全部相对论惯性的起因。一个明确给定的波动场模型可以出现类似的能量和惯性关系，但该模型已经具有相应的洛伦兹对称性。
- **开放的猜想：** 这些对应关系是否提示真实时空具有某种连续介质本性，仍然没有得到验证。该联想与真实世界没有进一步联系，也是完全可能的结果。

这是一份 **exploratory note / 思考笔记**，尚未经过独立的物理学同行评审，也未提交 arXiv。数值检查只核对指定模型中的数学关系，不构成物理实验或对时空猜想的验证。

建议先读论文第 3.3 节（延迟势的前后对称性）、第 7.2 节和附录 C（孤子模型），再读第 7.3 节（明确标记的猜想）。

## English summary

This note explores an intuitive association between the sound barrier and relativistic inertia. It follows the analogy through explicitly stated classical models and then considers a broader conjecture about effective continua and spacetime.

The retarded-field calculations do not establish gravitational pile-up as the origin of the light-speed limit. An acoustic metric illustrates the geometric description of wave propagation, and a sine-Gordon kink provides a concrete field-theoretic example with

$$
E = \gamma_{\ast} E_0, \qquad
P = \gamma_{\ast} \frac{E_0}{c_{\ast}^2}v, \qquad
\frac{\mathrm dP}{\mathrm dv} = \frac{E_0}{c_{\ast}^2}\gamma_{\ast}^3.
$$

These are properties of specified models, not evidence that the physical vacuum is a fluid. The manuscript separates its calculations from its speculative interpretation and identifies what additional dynamics and observational distinctions would be needed.

The author is a computer science student pursuing this idea in his spare time. OpenAI's GPT assisted substantially with mathematical development, literature lookup, organization, and drafting. The manuscript includes an explicit author and AI-use statement. This repository preserves an exploratory note and invites corrections; it does not report a validated new physical theory or new experimental results.

## Files

| File | Purpose |
| --- | --- |
| [Paper PDF](output/pdf/wave_speed_inertia.pdf) | English note with 2 figures, 3 appendices, and 12 references |
| [manuscript.tex](manuscript.tex) | Standalone LaTeX source; figures and bibliography are embedded |
| [verify_calculations.py](verify_calculations.py) | Deterministic consistency checks of the stated models |
| [verification_results.txt](verification_results.txt) | Output from the recorded verification run |
| [requirements.txt](requirements.txt) | Python dependency for reproducing those checks |

## Reproduce the calculations

Using Python 3 with NumPy:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verify_calculations.py
```

The script checks retarded-time geometry, the relativistic momentum derivative, the acoustic metric density, the angular transverse-traceless projection, the kink field equation and integrated energy–momentum, the quadrupole power, and the numerical estimates. It does not solve the nonlinear Einstein equations or test the continuum conjecture against observations.

## Build the PDF

With a LaTeX distribution containing the standard packages used by the source:

```sh
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
```

The reading PDF in this repository was exported with Tectonic. No external images, bibliography files, or shell escape are needed.

## Feedback

Corrections, relevant references, and suggestions for concrete tests are welcome through GitHub Issues. When reporting a mathematical issue, please identify the section or equation and the assumptions involved. The purpose of this repository is to keep a tentative idea accessible and open to revision.
