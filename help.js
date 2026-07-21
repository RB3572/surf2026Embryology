/* Shared "?" help system: a corner icon on each analysis opens a modal with an
 * accessible, figure-illustrated explanation. Pages add <button class="help-btn"
 * data-help="KEY">?</button>; the content lives in HELP[KEY]. Written in the plain-analogy
 * voice of the SURF deck (height-on-the-wall clock, the World Cup problem, campus buildings). */
(() => {
  // small inline SVG figures (fixed colors — the modal is always light)
  const BLUE = "#2563eb", RED = "#dc2626", ORANGE = "#e0752f", PURPLE = "#5c4d9e",
        GREEN = "#16a34a", INK = "#1f2937", GREY = "#94a3b8", PINK = "#db2777";

  const figPlanes = `<svg class="help-fig" viewBox="0 0 320 170" xmlns="http://www.w3.org/2000/svg">
    <circle cx="160" cy="85" r="72" fill="#f1f5fb" stroke="#dce3ee"/>
    <line x1="160" y1="20" x2="160" y2="150" stroke="${GREY}" stroke-width="1.5" stroke-dasharray="4 4"/>
    <line x1="96" y1="150" x2="224" y2="20" stroke="${ORANGE}" stroke-width="3"/>
    <circle cx="128" cy="60" r="8" fill="${BLUE}"/><circle cx="150" cy="45" r="8" fill="${BLUE}"/>
    <circle cx="185" cy="70" r="8" fill="${BLUE}"/><circle cx="200" cy="118" r="8" fill="${RED}"/>
    <circle cx="160" cy="85" r="5" fill="${INK}"/>
    <text x="235" y="30" font-size="11" fill="${ORANGE}">best line</text>
    <text x="235" y="128" font-size="11" fill="${INK}">3 vs 1</text>
  </svg>`;
  const figClock = `<svg class="help-fig" viewBox="0 0 320 160" xmlns="http://www.w3.org/2000/svg">
    <line x1="40" y1="18" x2="40" y2="150" stroke="${INK}" stroke-width="2"/>
    ${[0,1,2,3,4].map(i=>`<line x1="34" y1="${30+i*28}" x2="46" y2="${30+i*28}" stroke="${GREY}" stroke-width="1.5"/>`).join("")}
    <circle cx="120" cy="120" r="9" fill="${PURPLE}"/><circle cx="150" cy="118" r="9" fill="${PURPLE}"/>
    <text x="112" y="145" font-size="10" fill="${GREY}">early · far apart</text>
    <circle cx="235" cy="55" r="9" fill="${PURPLE}"/><circle cx="252" cy="55" r="9" fill="${PURPLE}"/>
    <text x="222" y="38" font-size="10" fill="${GREY}">late · close</text>
    <path d="M170 112 Q205 90 232 66" stroke="${GREEN}" stroke-width="2" fill="none" marker-end="url(#ar)"/>
    <defs><marker id="ar" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto"><path d="M0 0 L6 3 L0 6 z" fill="${GREEN}"/></marker></defs>
  </svg>`;
  const figSeg = `<svg class="help-fig" viewBox="0 0 320 170" xmlns="http://www.w3.org/2000/svg">
    <circle cx="120" cy="90" r="66" fill="#f6f4fb" stroke="#e4dcf1"/>
    <circle cx="120" cy="90" r="30" fill="#efe9fa" stroke="#d9ccf0"/>
    <circle cx="175" cy="52" r="15" fill="#fdeede" stroke="#f2d7b3"/>
    ${Array.from({length:22}).map(()=>{const a=Math.random()*6.28,r=8+Math.random()*20;return `<circle cx="${120+r*Math.cos(a)}" cy="${90+r*Math.sin(a)}" r="2.6" fill="${PURPLE}"/>`;}).join("")}
    ${Array.from({length:9}).map(()=>{const a=Math.random()*6.28,r=Math.random()*11;return `<circle cx="${175+r*Math.cos(a)}" cy="${52+r*Math.sin(a)}" r="2.6" fill="${ORANGE}"/>`;}).join("")}
    <text x="96" y="94" font-size="10" fill="#6d5aa8">pronuclei</text>
    <text x="150" y="30" font-size="10" fill="${ORANGE}">polar body</text>
    <text x="70" y="165" font-size="11" fill="${INK}">some segments are packed, others sparse</text>
  </svg>`;
  const figAxes = `<svg class="help-fig" viewBox="0 0 320 170" xmlns="http://www.w3.org/2000/svg">
    <circle cx="160" cy="90" r="66" fill="#fdeef3" stroke="#f6d3e0"/>
    <circle cx="160" cy="24" r="7" fill="${INK}"/><text x="172" y="22" font-size="10" fill="${INK}">polar body (animal pole)</text>
    <circle cx="215" cy="128" r="6" fill="${GREEN}"/><text x="150" y="150" font-size="10" fill="${GREEN}">sperm entry</text>
    <line x1="96" y1="70" x2="224" y2="110" stroke="${PINK}" stroke-width="3" stroke-dasharray="6 4"/>
    <text x="228" y="112" font-size="10" fill="${PINK}">1st cleavage?</text>
  </svg>`;
  const figPCA = `<svg class="help-fig" viewBox="0 0 320 170" xmlns="http://www.w3.org/2000/svg">
    <circle cx="150" cy="90" r="68" fill="#eefaf1" stroke="#cfeedd"/>
    ${Array.from({length:26}).map((_,i)=>{const a=-0.6+ (Math.random()-.5)*0.7, r=10+Math.random()*46;return `<circle cx="${150+r*Math.cos(a)}" cy="${90+r*Math.sin(a)}" r="2.6" fill="${BLUE}" opacity="0.8"/>`;}).join("")}
    <line x1="120" y1="108" x2="205" y2="66" stroke="${INK}" stroke-width="2.5" marker-end="url(#ar2)"/>
    <circle cx="212" cy="62" r="6" fill="${RED}"/><text x="220" y="60" font-size="10" fill="${RED}">sperm</text>
    <text x="70" y="120" font-size="10" fill="${BLUE}">gene cloud leans toward the sperm</text>
    <defs><marker id="ar2" markerWidth="9" markerHeight="9" refX="6" refY="3" orient="auto"><path d="M0 0 L6 3 L0 6 z" fill="${INK}"/></marker></defs>
  </svg>`;
  const figAB = `<svg class="help-fig" viewBox="0 0 320 160" xmlns="http://www.w3.org/2000/svg">
    <circle cx="110" cy="80" r="52" fill="#fbede2" stroke="#f0cfa9"/><text x="102" y="85" font-size="20" font-weight="700" fill="${ORANGE}">A</text>
    <circle cx="210" cy="80" r="52" fill="#ece9f6" stroke="#cfc6ea"/><text x="202" y="85" font-size="20" font-weight="700" fill="${PURPLE}">B</text>
    <text x="86" y="150" font-size="10" fill="${GREY}">sperm entry →</text>
    <text x="176" y="150" font-size="10" fill="${GREY}">← more transcripts?</text>
    <text x="150" y="20" font-size="11" fill="${INK}">which one is α?</text>
  </svg>`;

  const HELP = {
    "division-planes": { eyebrow: "Division Planes · the method", title: "The World Cup problem",
      html: `<p class="lede">Four superstars stand on the pitch. Draw <b>one</b> straight line to split them into
        two teams — but the line must pass through a fixed dot in the middle. Which line puts the most players on a
        single side? You can't reason it out; you try every angle and keep the most lopsided split.</p>
        ${figPlanes}<div class="help-cap">Every line pivots about the same center dot; we keep the most unbalanced split.</div>
        <h3>What we actually measure</h3>
        <p>In the zygote the "players" are a gene's <b>transcripts</b> and the fixed dot is the cell's
        <b>center of mass</b>. A line becomes a <b>plane</b> that must contain the <b>polar-body axis</b>
        (center → polar body). We sweep <b>18 planes</b>, 10° apart, and for each one count how many of that gene's
        transcripts fall on side A vs side B.</p>
        <h3>How to read the 3-D view</h3>
        <ul><li><span class="tag">blue</span> / <span class="tag">red</span> dots — the selected gene's transcripts on side A / side B of the shown plane.</li>
        <li><span class="tag">green</span> dots — molecules outside the cytoplasm (polar body, pronuclei); not counted.</li>
        <li>The orange square is the current plane; "All planes" colors every plane by significance (dark = a strong, lopsided split).</li></ul>
        <h3>Why it matters</h3>
        <p>A gene whose transcripts are <b>consistently lopsided</b> across the best plane hints that the zygote is
        already <b>pre-patterned</b> — that molecular asymmetry, present before the first division, could mark where
        the two-cell blastomeres (α and β) will differ.</p>
        <div class="help-callout">The <b>p-value</b> asks: could a split this lopsided happen by chance if each
        molecule flipped a fair coin for its side? Low p = unlikely to be luck.</div>` },

    "division-crossembryo": { eyebrow: "Division Planes · bottom drawer", title: "Do all zygotes lean the same way?",
      html: `<p class="lede">One zygote leaning to one side could be luck. The real question is whether <b>many</b> zygotes
        lean the <b>same</b> way once you line them up.</p>
        <h3>The two panels</h3>
        <ul><li><b>Aligned cross-sections</b> — each zygote's outline, rotated so its own best plane is vertical and flipped so
        the higher-count side is on the right. Each is colored by how significant its split is (dark purple = significant).</li>
        <li><b>Per-side counts</b> — for the chosen gene, a stacked bar of left-vs-right transcript counts per zygote. A consistent
        lean across embryos is the signal.</li></ul>
        <h3>The panel caveat</h3>
        <p>The zygotes come from several MERFISH gene panels, so <b>no single gene is measured in all of them</b>. Picking an
        alignment gene shows only the zygotes that contain it — that's why the count changes with the gene.</p>` },

    "circularize": { eyebrow: "Division Planes · toggle", title: "Circularize the cell",
      html: `<p class="lede">Real cells are lumpy, and a lump can fake a lopsided split. So we <b>inflate</b> each cell
        until it's a smooth sphere — then ask whether the asymmetry <i>survives</i>.</p>
        <div class="help-cap">A shape control: normalize the cell's shape, keep the molecules' relative positions.</div>
        <h3>What the toggle does</h3>
        <p>For the cytoplasm (segment 1) only, every point is pushed onto a smooth sphere of the cell's <b>average
        radius</b>. Transcripts move <b>with</b> the tissue (they stay inside, they don't pile onto the shell). Then the
        <b>entire</b> analysis — center of mass, per-side counts, p-values, cross-section — is recomputed on the rounded cell.</p>
        <h3>Why</h3>
        <p>If a gene's split stays lopsided after the shape is normalized, the asymmetry is about <b>where the molecules are</b>,
        not about the cell being an odd shape. It's a shape-control for the pre-patterning question.</p>` },

    "pronuclei": { eyebrow: "Pronuclei Distance · the model", title: "A clock you can read from a ruler",
      html: `<p class="lede">Mark a child's height on the wall each year. The marks creep upward, and the <b>spacing</b> becomes
        a clock — you can read <i>age</i> from <i>height</i>. We do the same with the zygote, reading developmental <b>time</b>
        from a physical <b>distance</b>.</p>
        ${figClock}<div class="help-cap">After fertilization the two pronuclei form near the surface and migrate together; the gap shrinks as the zygote ages.</div>
        <h3>The measurement</h3>
        <p>The two pronuclei (maternal + paternal) start far apart and converge just before the first division. So the
        <b>minimum distance between them</b> is a stand-in for developmental time — <b>smaller distance = later</b>. We plot each
        zygote's transcript count against that distance and fit a curve.</p>
        <h3>How to read it</h3>
        <ul><li>Each dot is one zygote; x = pronuclei distance (µm), y = transcript count.</li>
        <li>The <b>Regression model</b> menu fits different shapes (linear, exponential decay of maternal mRNA, a ZGA-style
        sigmoid, count-GLMs, …) — a way to ask <i>how</i> count changes with time, not just whether it does.</li>
        <li><b>Region</b> restricts the count to one part of the embryo — the whole embryo, segment 1 (cytoplasm), the two
        <b>pronuclei</b> (auto-detected, the same segments used for the distance), or the <b>polar bodies</b> (every labelled
        segment that is not segment 1 or a pronucleus).</li>
        <li><b>Count axis</b> plots the raw count, the count ÷ the zygote's total transcripts (a fraction), or ÷ the region's
        volume (density per µm³ — controls for a bigger region simply holding more molecules).</li>
        <li><b>Pseudotime</b> flips the x-axis to (max distance − distance) so larger = later; <b>Flip x / y</b> transposes any
        graph (the fit is unchanged).</li></ul>
        <h3>Gene sets</h3>
        <p>The bottom graph sums a <b>set</b> of genes' counts per zygote and plots that against distance. Build a set by hand or
        from a <b>preset</b>; presets <i>add</i> to the current set. The two highlighted presets are <b>live</b> — <b>Top 10
        ＋/−correlated</b> pull the current strongest positive / negative distance-correlates straight from the ranking. The rest
        are curated: <b>functional clusters</b> whose members individually track distance and share a role (Notch/Wnt/Hedgehog,
        Ras–MAPK, proteostasis, oocyte/pluripotency TFs, "rises toward first cleavage"), plus biology sets from the deck
        (maternally-deposited, ZGA markers, pronucleus-associated, …).
        By default a zygote is plotted if it contains <i>any</i> set gene (missing genes are simply left out of its sum); tick
        <b>require all genes</b> for a strictly comparable sum over an identical gene list — best for a set drawn from one MERFISH
        panel, since the panels are disjoint.</p>
        <div class="help-callout"><b>Total</b> transcript count is a weak clock (r ≈ 0.16). Individual genes — and curated gene
        sets — do far better; the right drawer ranks single genes by that correlation.</div>` },

    "pronuclei-genes": { eyebrow: "Pronuclei · gene ranking", title: "Which single gene is the best clock?",
      html: `<p class="lede">The whole transcriptome is a blurry clock. But a <b>single well-behaved gene</b> can be sharp — its
        count might track pronuclei distance almost perfectly.</p>
        <h3>What the ranking shows</h3>
        <p>For every gene present in enough zygotes, we measure the <b>Pearson correlation</b> between its transcript count and the
        pronuclei distance, and rank the strongest <b>positive</b> and <b>negative</b> relationships. Positive = more of the gene
        goes with a larger gap; negative = more of the gene goes with a smaller gap.</p>
        <p>Click a gene to see its own scatter above the total, and toggle its molecules as dots in the 3-D cell. Because panels
        are disjoint, each gene is measured in a different subset of zygotes (its <span class="tag">n</span>).</p>` },

    "pronuclei-stats": { eyebrow: "Pronuclei · the statistics", title: "What every number under a graph means",
      html: `<p class="lede">Each scatter (selected gene, all transcripts, and the gene set) reports the same five things, so you
        can judge at a glance whether a trend is real or could be luck.</p>
        <svg class="help-fig" viewBox="0 0 320 150" xmlns="http://www.w3.org/2000/svg">
          <line x1="30" y1="18" x2="30" y2="128" stroke="${INK}" stroke-width="1.5"/>
          <line x1="30" y1="128" x2="150" y2="128" stroke="${INK}" stroke-width="1.5"/>
          <text x="16" y="88" font-size="9" fill="${GREY}" transform="rotate(-90 16 88)">count</text>
          <text x="60" y="144" font-size="9" fill="${GREY}">distance</text>
          ${[[42,110],[60,96],[74,100],[92,78],[108,70],[128,52]].map(p=>`<circle cx="${p[0]}" cy="${p[1]}" r="4" fill="${BLUE}"/>`).join("")}
          <line x1="38" y1="112" x2="140" y2="52" stroke="${ORANGE}" stroke-width="2.5"/>
          <text x="150" y="30" font-size="10" fill="${ORANGE}">y = m·x + b</text>
          <text x="170" y="70" font-size="9" fill="${GREY}">real: points hug the line</text>
          <line x1="180" y1="88" x2="300" y2="88" stroke="#e5e7eb" stroke-width="1"/>
          <text x="170" y="104" font-size="9" fill="${GREY}">null: shuffle count↔distance</text>
          ${[[190,118],[210,100],[228,120],[246,102],[266,116],[288,104]].map(p=>`<circle cx="${p[0]}" cy="${p[1]}" r="3" fill="${GREY}"/>`).join("")}
        </svg>
        <div class="help-cap">The fit line is the model; the p-value asks how often a random re-pairing of the same points would look as trend-y.</div>
        <h3>The five numbers</h3>
        <ul>
        <li><b>Equation</b> — the fitted curve's formula. For the default <b>Linear</b> model it is <b>y = m·x + b</b>: the slope
        <b>m</b> is transcripts gained (or lost) per micron of pronuclei distance, and <b>b</b> is the intercept. Other models show
        their own form (e.g. <span class="tag">y = a·e^(k·x)</span> for exponential).</li>
        <li><b>R²</b> — the share of the variation the fit explains (1 = perfect, 0 = no better than a flat line). Because some
        models are fitted in a transformed space, R² is reported on that model's <b>natural scale</b> — labelled <span class="tag">R²</span>,
        <span class="tag">R² (log)</span>, <span class="tag">R² (logit)</span> or <span class="tag">pseudo-R²</span> — so a
        log-fit's R² isn't compared unfairly against a raw one.</li>
        <li><b>Pearson r</b> — the straight-line correlation of distance and count, from <b>−1</b> to <b>+1</b>. The sign is the
        direction (＋ = count rises with distance); |r| is the strength. r² equals the linear R².</li>
        <li><b>p-value</b> — the chance of seeing a correlation this strong <b>if distance and count were actually unrelated</b>
        (two-sided). Small p = unlikely to be a fluke. We compute it two independent ways that agree:</li>
        </ul>
        <h3>How the p-value is calculated</h3>
        <p><b>1 · Exact test.</b> Turn r into a t-statistic, <b>t = r·√((n−2)/(1−r²))</b>, which follows a Student-t distribution
        with <b>n−2</b> degrees of freedom when there is no real correlation; the two-sided tail area is the p-value (via the
        incomplete-beta function). This is the classic significance test for a correlation.</p>
        <p><b>2 · Permutation null (the honest cross-check).</b> Keep the same numbers but <b>shuffle which count is paired with
        which distance</b> 2000 times — this is a world where, by construction, there is <i>no</i> relationship. Each shuffle gets
        its own r; the p-value is how often that random |r| is at least as big as the one we actually saw. Both methods land in the
        same place (shown together in the tooltip on the <span class="tag">p</span>).</p>
        <p>Shorthand stars: <b>***</b> p&lt;0.001, <b>**</b> p&lt;0.01, <b>*</b> p&lt;0.05, <b>ns</b> = not significant.</p>
        <div class="help-callout"><b>Watch two traps.</b> (1) <b>n</b> is small and varies per gene (disjoint panels), so one
        outlier can swing r — trust <b>*</b>/<b>**</b> more when n is larger. (2) The p-value tests a <i>straight-line</i>
        association; a curved model can have a high R² yet a near-zero Pearson r (a U-shape) — so read the R², the r, <i>and</i> the
        shape together, not any one alone.</div>` },

    "extpt": { eyebrow: "Extended Pseudotime · the idea", title: "One clock across three stages",
      html: `<p class="lede">The pronuclei give zygotes a ruler for developmental time. This project extends that ruler
        past the first division — laying <b>zygote → early-2-cell → late-2-cell</b> along a single pseudotime axis so you
        can watch a gene's transcripts rise and fall across the whole maternal-to-zygotic transition.</p>
        <h3>The axis</h3>
        <p>Three stage blocks in developmental order (zygote 0–1, early-2C 1–2, late-2C 2–3), coloured and banded. Within
        the zygote block, embryos are ordered by the <b>pronuclei-distance clock</b> (larger gap = earlier). Points are
        coloured by stage; the fitted model runs across all three.</p>
        <h3>Is there a clock for the 2-cell embryos?</h3>
        <p>Not a geometric one — once the cell divides there are no pronuclei to measure. So within the early- and
        late-2-cell blocks we order embryos by <b>total transcript count</b> as a maturation proxy (the transcriptome
        grows as zygotic genome activation ramps up). Treat this as a <i>soft</i> ordering: the honest signal is the
        <b>cross-stage trajectory</b> (how the three bands step up or down), not the fine ordering inside a 2-cell block.
        Because that block is ordered by total count, the "All transcripts" graph is diagonal within it by construction.</p>
        <h3>Same machinery as Pronuclei</h3>
        <p>Per-gene, all-transcripts and gene-set graphs; the full regression-model menu; complete statistics (R²,
        equation, Pearson r, and a null-tested p-value — see the <b>?</b> by the model menu); tabs, corner-resize and
        high-res download. The right drawer ranks genes by their correlation with pseudotime across all three stages
        (405 genes appear in every stage, so trajectories are data-rich).</p>` },

    "segments": { eyebrow: "Segment Enrichment · the idea", title: "People on campus",
      html: `<p class="lede">At any moment some campus buildings are packed and others nearly empty. If you know which building
        someone is in, you often know what they're doing. Molecules are the same: <b>where</b> a transcript sits often predicts
        <b>what</b> it's for.</p>
        ${figSeg}<div class="help-cap">Every transcript is assigned to a segment — cell body, polar body, pronuclei — then we ask which genes are over-represented where.</div>
        <h3>What "enriched" means</h3>
        <p>For each segment we compute a gene's <b>density fold-change</b>: how concentrated it is in that segment versus the
        cell-wide average, correcting for the segment's volume. A gene wholly inside one small segment scores very high. (You can
        switch to a simpler <b>fraction</b> mode: what share of the gene's molecules land in the segment.)</p>
        <h3>Why it matters</h3>
        <p>Subcellular localization is real and pervasive — in fly embryos, 71% of expressed genes had distinct patterns
        (Lécuyer et al., Cell 2007). Genes that concentrate in the pronuclei or polar body point to the machinery running those
        structures.</p>` },

    "axes": { eyebrow: "Fertilization Geometry · the question", title: "Do two landmarks predict the first cut?",
      html: `<p class="lede">Give someone two fixed landmarks on a balloon — where you poked it, and its north pole — and ask them
        to guess where it will first split. That's the pre-patterning question for the embryo.</p>
        ${figAxes}<div class="help-cap">Sperm entry (the poke) + the polar body (the animal pole) — do they predict the first cleavage plane?</div>
        <h3>The landmarks</h3>
        <ul><li><b>Sperm entry</b> — the GFP-labeled midpiece marks where the sperm fused. A proposed organizer of the first axis.</li>
        <li><b>Polar body</b> — the extruded body sits at the <b>animal pole</b>; the opposite side is vegetal.</li>
        <li><b>First cleavage plane</b> — read directly in two-cell embryos (the interface between blastomeres); in zygotes we use
        the pronuclear axis as the predicted normal.</li></ul>
        <h3>How it's tested</h3>
        <p>We measure the angle between these landmarks and the cleavage plane and compare it to an <b>exact random-orientation
        null</b> — the pattern you'd see if orientation were pure chance. A <b>shape control</b> also checks the signal isn't just
        the cell being elongated. Modest n, hypotheses fixed in advance.</p>` },

    "sperm-pca": { eyebrow: "Sperm Prediction · the idea", title: "Reading the sperm's wake",
      html: `<p class="lede">A boat leaves a wake that points back to where it went. Some genes' transcript clouds seem to
        <b>lean</b> toward the sperm — if enough of them agree, their combined "arrow" can <b>predict</b> where the sperm is.</p>
        ${figPCA}<div class="help-cap">Each gene's cloud has a long axis (PCA) and a center of mass; consistent genes point at the sperm.</div>
        <h3>What we compute</h3>
        <p>For each gene we take its transcript cloud's <b>principal axis</b> (PCA) and its <b>center of mass</b>, and measure how
        well they align with the direction to the manually-located sperm. Genes are <b>ranked by consistency</b> across embryos;
        the top genes' directions are averaged into a predicted sperm location.</p>
        <h3>How to read it</h3>
        <p>A <b>randomized null</b> (shuffled sperm positions) shows what accuracy chance alone would give — the real prediction
        has to beat it to mean anything.</p>` },

    "alphabeta": { eyebrow: "Sperm α/β · the chart", title: "Which twin is α?",
      html: `<p class="lede">The two-cell embryo has two blastomeres that look almost identical. Several independent "tests"
        each nominate one of them as <b>alpha</b> — but a test only says "this one," not whether "this one" is the same twin
        another test meant. The chart lines the tests up so you can see if they <b>agree</b>.</p>
        ${figAB}<div class="help-cap">A = original blastomere labels 1+3, B = labels 2+4. Each row is a labeling method.</div>
        <h3>The methods</h3>
        <p><b>Exact:</b></p>
        <ul><li><b>Sperm entry</b> — the blastomere the sperm entered is called alpha.</li>
        <li><b>Higher total transcript</b> — the blastomere with more transcripts is alpha (matches Harry's grids 20/20).</li></ul>
        <p>Rows marked <b style="color:#db2777">≈</b> are <b>best-guess reconstructions</b> of Harry's methods — his deck
        names each and shows results, but not the exact parameters, so these approximate the likely implementation:</p>
        <ul><li><b>Expression-axis PCA</b> — project each blastomere's gene profile onto the dominant axis of expression variation.</li>
        <li><b>Decreased / Increased panel</b> — genes split into <i>maternal</i> (fade from zygote to 2-cell) vs <i>zygotic</i> (rise).</li>
        <li><b>… : ratio-sum</b> — add up each gene's volume-normalized A-vs-B log-ratio over that panel; the sign picks alpha.</li>
        <li><b>… : PCA</b> — the same per-gene asymmetries, but weighted toward the dominant shared pattern instead of counted evenly.</li>
        <li><b>Exhaustive</b> — score each gene by how far its split departs from an even one, then let the most decisive genes vote
        (<i>unfiltered</i> = all genes; <i>mean count ≥ 20</i> = only well-expressed genes).</li></ul>
        <h3>Why the flips</h3>
        <p>"Alpha vs beta" is just a name, and each method may pick the name in its own arbitrary direction. So each row can be
        <b>flipped</b> (swap A↔B). <b>Auto-align</b> flips rows to make the columns as uniform as possible; a column where every
        method agrees is a <b>confident consensus</b> call for that embryo. <b>Flip all</b> swaps the global sense; the
        <b>consensus</b> row summarizes each column.</p>
        <div class="help-callout">The goal of all of this: identify the division plane in the zygote and tell the two-cell
        blastomeres apart — the same asymmetry, one stage later.</div>` },
  };

  // ---------- modal machinery ----------
  let overlay;
  function build() {
    overlay = document.createElement("div");
    overlay.className = "help-overlay";
    overlay.innerHTML = `<div class="help-modal" role="dialog" aria-modal="true">
      <div class="help-head"><div><div class="help-eyebrow" id="help-eyebrow"></div>
        <div class="help-title" id="help-title"></div></div>
        <button class="help-x" id="help-x" aria-label="Close">✕</button></div>
      <div class="help-body" id="help-body"></div></div>`;
    document.body.appendChild(overlay);
    overlay.addEventListener("click", (e) => { if (e.target === overlay) close(); });
    overlay.querySelector("#help-x").addEventListener("click", close);
    document.addEventListener("keydown", (e) => { if (e.key === "Escape") close(); });
  }
  function open(key) {
    const c = HELP[key]; if (!c) return;
    if (!overlay) build();
    overlay.querySelector("#help-eyebrow").textContent = c.eyebrow || "";
    overlay.querySelector("#help-title").textContent = c.title || "";
    overlay.querySelector("#help-body").innerHTML = c.html || "";
    overlay.querySelector(".help-modal").scrollTop = 0;
    overlay.classList.add("open");
  }
  function close() { if (overlay) overlay.classList.remove("open"); }

  // delegate clicks from any [data-help] button (works for buttons added later too)
  document.addEventListener("click", (e) => {
    const b = e.target.closest("[data-help]");
    if (b) { e.preventDefault(); e.stopPropagation(); open(b.getAttribute("data-help")); }
  });
  window.Help = { open, close, has: (k) => k in HELP };
})();
