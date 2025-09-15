
# Preface0

For many students, their first real encounter with spectral analysis often unfolds like this:

One day, they notice an intriguing phenomenon in a time-domain signal and eagerly share their discovery with an advisor or senior colleague. The response is calmly delivered: "You should try some Fourier or wavelet analysis."

Returning to their desk, they dig out an old calculus textbook, flip through several pages, and quickly realize it won’t help. Next comes the trusted solution: a swift Google search for “Fourier Analysis tutorial". Eventually, they stumble upon a highly recommended *Digital Signal Processing* textbook with an impressive 4.3/5.0 book rating.  After a marathon weekend, they manage to read through the 50-plus pages of Chapter 1, *Signals and Systems*. However, by Chapter 2, Linear Time-Invariant Systems, fatigue sets in—only to realize that the actual Fourier series material is still more than 120 pages away.

At this juncture, most students pragmatically pivot to google *"Fourier Analysis by xxx"*  and get an answer with some unfamiliar jargon from *StackOverflow*, grabbing a ready-made code snippet to forge ahead.

Yet, a few determined souls persist—spending days gathering materials, watching lectures online, coding, and compiling a detailed report. Proudly, they present their hard work to their advisor, only to be met with the classic understated response: “Why so little progress this week?”

<div STYLE="page-break-after: always;"></div>

```mermaid
graph LR

subgraph L0["Level 0: $$[{f_s}/{4}, {f_s}/{2})$$"]
A --> B --Downsampling<br>↓2--> C
A --> D --↓2--> E
end

subgraph L1["Level 1: $$[{f_s}/{8}, {f_s}/{4})$$"]
E --> F --↓2--> G
E --> H --↓2--> I
end

subgraph L2["Level 2: $$[{f_s}/{16}, {f_s}/{8})$$"]
I --> J --↓2--> K
I --> L --↓2--> M
end

A@{ shape: lean-r, label: "$$x_0[n] = x[n]$$" }
B["$$x_0[n]*h[n]$$"]
C@{ shape: lean-l, label: "Level 0 Details<br>(cD0)" }
D["$$x_0[n]*g[n]$$"]
E@{ shape: lean-r, label: "$$x_1[n]$$" }
F["$$x_1[n]*h[n]$$"]
G@{ shape: lean-l, label: "Level 1 Details<br>(cD1)" }
H["$$x_1[n]*g[n]$$"]
I@{ shape: lean-r, label: "$$x_2[n]$$" }
J["$$x_2[n]*h[n]$$"]
K@{ shape: lean-l, label: "Level 2 Details<br>(cD2)" }
L["$$x_2[n]*g[n]$$"]
M@{ shape: lean-l, label: "Level 2 Approximate<br>(cA)" }

B@{shape: stadium}
D@{shape: stadium}
F@{shape: stadium}
H@{shape: stadium}
J@{shape: stadium}
L@{shape: stadium}

classDef lowpass  fill:#C5E1A5,stroke:#558B2F,color:#000,shape:stadium
classDef hipass   fill:#EF9A9A,stroke:#B71C1C,color:#000,shape:stadium
classDef coeff    fill:#BBDEFB,stroke:#1E88E5,color:#000,stroke-width:2px,shape:stadium

class B,F,J lowpass
class D,H,L hipass
class A,C,E,G,I,K,M coeff
```

```mermaid
graph LR

A(["$$x[n]$$"]) -->|Window|B@{ shape: rect, label: "$$x_{windowed}[n]$$" }

A-->|Padding|C@{shape:rect, label: "$$x_{padded}[n]$$"}

C -->|Window|B

B -->|Fourier
Transform|D@{shape: rect, label: "$$X[k]$$"}
```

