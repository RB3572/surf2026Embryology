/* gene-info.js — global gene hover tooltips.
 *
 * A curated map of concise, EARLY-EMBRYO-context descriptions (zygote → 2-cell → morula),
 * plus a self-running auto-tagger that adds a native `title` (hover) tooltip to any gene
 * shown anywhere on the site — dropdown options, ranking rows, gene-name spans, readouts —
 * and keeps tagging as pages render gene lists dynamically (MutationObserver).
 *
 * Coverage is deliberately partial: only genes with well-established roles in the early
 * mouse embryo are annotated; unfamiliar genes simply get no tooltip (never a guessed one).
 * To extend, add an entry to DESC — the key is the exact gene symbol as displayed.
 */
window.GeneInfo = (function () {
  "use strict";

  // Each blurb: core function · early-embryo role · what asymmetry between blastomeres may mean.
  const DESC = {
    // ── maternal-effect / subcortical maternal complex (SCMC) / oocyte store ──
    Nlrp5: "Subcortical maternal complex (SCMC) core (Mater). Maternally deposited; essential for progression past the 2-cell stage. High in the zygote, cleared as the zygotic program takes over.",
    Padi6: "SCMC component; builds the oocyte cytoplasmic lattices. Maternally deposited; required for ZGA and cleavage — loss arrests at 2-cell.",
    Nlrp2: "SCMC-associated maternal-effect gene. Maternally deposited; its loss causes early cleavage arrest and imprinting defects.",
    Nlrp9c: "Maternal NLRP of the oocyte store (SCMC-related).",
    Tle6: "Subcortical maternal complex member. Maternal; needed for the first cleavage divisions.",
    Zar1: "Maternal-effect regulator of the oocyte-to-embryo transition. Maternally deposited; loss blocks development at the 1–2-cell stage.",
    Mos: "Oocyte Mos–MEK–MAPK kinase. Maternal; maintains meiotic (CSF) arrest and is cleared after fertilization.",
    Fbxo43: "Emi2 — APC/C inhibitor. Maternal; holds the metaphase-II (CSF) arrest until fertilization triggers its destruction.",
    Dnmt1: "Maintenance DNA methyltransferase (oocyte Dnmt1o). Maternal; preserves genomic imprints through cleavage divisions.",
    Uhrf1: "Targets Dnmt1 to hemimethylated DNA. Maternal contribution maintains methylation in early cleavage.",
    Zp1: "Zona pellucida glycoprotein. Oocyte-made (maternal); part of the sperm-binding coat; transcripts decline after fertilization.",
    Zp2: "Zona pellucida glycoprotein and sperm receptor. Oocyte-made (maternal); cleaved at fertilization to block polyspermy; message depletes afterward.",
    Zp3: "Zona pellucida glycoprotein / primary sperm receptor. Oocyte-made (maternal); declines after fertilization.",
    Btbd18: "piRNA-pathway factor of germ cells (maternal store).",
    Lhx8: "Oocyte-specific transcription factor of early folliculogenesis (maternal).",
    Fmn2: "Actin-nucleating formin. Maternal; drives the cytoplasmic actin mesh and off-centre spindle positioning in oocyte/zygote — an asymmetric-division factor.",
    Tcl1b2: "TCL1-family PI3K/AKT co-activator; maternal, oocyte-enriched.",
    Tcl1b3: "TCL1-family PI3K/AKT co-activator; maternal, oocyte-enriched.",
    Tcl1b4: "TCL1-family PI3K/AKT co-activator; maternal, oocyte-enriched.",

    // ── zygotic genome activation (ZGA) / totipotency ──
    Duxf1: "Dux family — a master ZGA activator. Transiently switches on the minor-ZGA / MERVL 2-cell program, then is silenced.",
    Duxf3: "Dux family — a master ZGA activator. Transiently induces the 2-cell / MERVL program at genome activation.",
    Zscan4a: "2-cell / '2C-like' totipotency marker. Transient at ZGA; supports telomere elongation and genome stability.",
    Zscan4b: "2-cell / '2C-like' totipotency marker; transient at ZGA (telomere/genome stability).",
    Zscan4d: "2-cell / '2C-like' totipotency marker; transient at ZGA (telomere/genome stability).",
    Zscan4e: "2-cell / '2C-like' totipotency marker; transient at ZGA (telomere/genome stability).",
    Zscan4f: "2-cell / '2C-like' totipotency marker; transient at ZGA (telomere/genome stability).",
    Obox1: "Oocyte-specific homeobox transcription factor; activates zygotic genome activation genes.",
    Obox2: "Oocyte-specific homeobox transcription factor; a driver of zygotic genome activation.",
    Obox3: "Oocyte-specific homeobox transcription factor; activates ZGA genes.",
    Obox8: "Oocyte-specific homeobox transcription factor; drives zygotic genome activation.",
    "MuERV-L": "MERVL endogenous retrovirus — the definitive 2-cell / totipotency marker; transiently expressed at ZGA, then silenced.",
    L1td1: "LINE-1-associated RNA-binding protein linked to the totipotent 2-cell program and pluripotency.",
    Zfp352: "2-cell-specific ZGA gene and Dux target.",
    Nr5a2: "Nuclear receptor (LRH-1) that helps trigger zygotic genome activation.",
    Gadd45a: "Couples to active DNA demethylation and the stress response around ZGA.",
    Kdm4dl: "KDM4/JmjC-family histone-demethylase relative expressed around ZGA.",
    Eif1ad12: "Translation-associated factor in the ZGA-linked gene cluster.",
    Trib3: "Pseudokinase; among genes induced at the 2-cell transition.",
    LincGET: "2-cell nuclear lncRNA that recruits CARM1. Its ASYMMETRY between blastomeres biases the first lineage decision toward the inner (pluripotent) cell fate.",
    Carm1: "Arginine methyltransferase (H3R26me2). Higher CARM1 in a 2–4-cell blastomere biases its progeny toward the pluripotent inner-cell-mass fate — a key early asymmetry.",
    Pramef8: "Preferentially-expressed-antigen family; 2-cell/ZGA-associated.",

    // ── pluripotency, ICM vs trophectoderm, endoderm ──
    Sox2: "Pluripotency transcription factor; the first TF to mark inner-cell-mass fate. Early Sox2 heterogeneity among blastomeres foreshadows ICM vs trophectoderm.",
    Klf2: "Naive-pluripotency transcription factor of the inner cell mass.",
    Klf4: "Naive-pluripotency transcription factor (inner cell mass / epiblast).",
    Esrrb: "Naive-pluripotency transcription factor (ICM / epiblast).",
    Tfcp2l1: "Naive-pluripotency transcription factor.",
    Foxd3: "Pluripotency / inner-cell-mass transcription factor.",
    Lin28a: "RNA-binding pluripotency factor; blocks let-7 microRNA maturation.",
    Myc: "Proliferation / biosynthesis driver; part of the CARM1–Myc axis that biases blastomeres toward the inner cell mass.",
    Yap1: "Hippo-pathway effector. Nuclear YAP in outer cells switches on trophectoderm genes; its inside–outside asymmetry drives the trophectoderm-vs-ICM decision.",
    Tead2: "TEAD transcription factor; partners YAP to specify trophectoderm.",
    Lats1: "Hippo kinase, active in inner cells; excludes YAP from the nucleus, favouring inner-cell-mass fate.",
    Gata3: "Trophectoderm transcription factor downstream of TEAD4/YAP.",
    Gata4: "Primitive-endoderm transcription factor (later inner-cell-mass segregation).",
    Sox7: "Primitive-endoderm marker.",
    Hand1: "Trophoblast-giant-cell transcription factor.",
    Id1: "BMP-target HLH factor that biases cells away from the pluripotent program.",
    Id4: "BMP-target HLH factor; differentiation-associated.",
    Cebpa: "Differentiation-associated transcription factor.",
    Cebpb: "Differentiation-associated transcription factor.",
    Elf3: "Epithelial ETS transcription factor associated with trophectoderm.",
    Klf6: "Differentiation-linked KLF transcription factor.",
    Tfap2a: "AP-2 transcription factor; epithelial / trophectoderm programs.",
    Tcf7: "Wnt/β-catenin transcription factor (TCF1).",

    // ── polarity, junctions, compaction, cleavage cytoskeleton ──
    Pard3: "Apical polarity scaffold; organizes tight junctions. From the 8-cell stage it establishes apical–basal polarity and inside–outside fate. Earlier asymmetry may indicate stronger or earlier polarity in some blastomeres.",
    Pard6a: "Par-complex polarity adaptor (with aPKC); apical-domain formation at polarization.",
    Prkci: "Atypical PKC (aPKC) — the polarity kinase of the Par complex; needed for apical-domain formation and trophectoderm specification.",
    Prkce: "Protein kinase C-epsilon; PKC signalling in the early embryo.",
    Ctnnb1: "β-catenin — E-cadherin adhesion at compaction plus Wnt signalling. Adhesion differences influence cell positioning and sorting.",
    Cdc42: "Rho-family GTPase that builds the apical domain during polarization.",
    RhoA: "Rho GTPase controlling actomyosin contractility for cleavage and compaction.",
    Rhou: "Atypical Rho GTPase (Wrch1) acting on junctions and actin.",
    Myh9: "Non-muscle myosin IIA; cortical tension for compaction and apical-domain positioning.",
    Fmnl2: "Formin-family actin regulator (motility and junctions).",
    Arpc1b: "Arp2/3-complex subunit; branched cortical actin networks.",
    Numb: "Cell-fate determinant that can segregate asymmetrically at division and inhibits Notch.",
    Krt8: "Keratin-8 intermediate filament; an early trophectoderm epithelial marker (TROMA-1).",
    Krt18: "Keratin-18 intermediate filament; trophectoderm epithelial marker.",
    Dsg2: "Desmoglein-2 desmosomal junction protein; trophectoderm epithelialization.",
    Pkp2: "Plakophilin-2 desmosomal protein; epithelial junctions.",
    Gja4: "Connexin gap-junction protein; intercellular coupling during compaction.",
    Gjb5: "Connexin gap-junction protein; intercellular coupling in the early embryo.",
    Itgb3: "Integrin β3; cell–matrix adhesion.",
    Sdc4: "Syndecan-4 adhesion/co-receptor.",
    Specc1l: "Cytoskeletal cross-linker supporting adhesion and junction integrity.",

    // ── signalling: Wnt · Notch · Nodal/TGF-β/BMP · FGF · Hedgehog ──
    Fzd2: "Frizzled Wnt receptor.",
    Fzd4: "Frizzled Wnt receptor.",
    Fzd5: "Frizzled Wnt receptor.",
    Axin2: "Wnt/β-catenin target and negative-feedback scaffold — a readout of pathway activity.",
    Wif1: "Secreted Wnt inhibitor.",
    Wnt7b: "Wnt ligand.",
    Nrarp: "Notch target that also tunes Wnt — a signalling-feedback node.",
    Zbed3: "Axin-binding Wnt-pathway activator.",
    Notch2: "Notch receptor; Notch signalling is active in the trophectoderm lineage.",
    Notch3: "Notch receptor.",
    Jag2: "Jagged Notch ligand.",
    Dtx2: "Deltex, a Notch-pathway regulator.",
    Hes1: "Notch-target bHLH repressor.",
    Nodal: "Nodal (TGF-β/Activin) ligand; establishes embryonic axes and supports pluripotency signalling.",
    Lefty2: "Nodal antagonist; shapes Nodal signalling gradients.",
    Smad1: "BMP-branch SMAD transcriptional effector.",
    Smad2: "Nodal/Activin-branch SMAD effector.",
    Smad7: "Inhibitory SMAD; negative feedback on TGF-β/BMP signalling.",
    Bmp5: "BMP ligand; promotes differentiation programs.",
    Bmp7: "BMP ligand; promotes differentiation / trophoblast programs.",
    Bmpr1b: "BMP type-I receptor.",
    Acvr2b: "Activin/Nodal type-II receptor.",
    Amhr2: "TGF-β-superfamily (AMH) type-II receptor.",
    Bambi: "BMP/TGF-β pseudoreceptor (decoy) that dampens signalling.",
    Inha: "Inhibin subunit; antagonizes activin signalling.",
    Fst: "Follistatin; binds and inhibits activin/BMP.",
    Fkbp1a: "FKBP12; keeps TGF-β-family receptors inactive until ligand binds.",
    Tgfb1: "TGF-β ligand.",
    Tgfb2: "TGF-β ligand.",
    Fgf4: "The key ICM-derived FGF; drives primitive-endoderm vs epiblast segregation in the inner cell mass.",
    Fgf8: "FGF ligand (patterning).",
    Egfr: "EGF receptor; growth-factor signalling.",
    Gli3: "Hedgehog-pathway transcription factor (repressor/activator).",
    Gpr161: "Ciliary GPCR that negatively regulates Hedgehog signalling.",
    Tulp3: "Ciliary trafficking factor; negative regulator of Hedgehog.",

    // ── chromatin / transcription / RNA regulation ──
    Setd2: "H3K36me3 methyltransferase; transcription-coupled chromatin mark.",
    Ep300: "p300 histone acetyltransferase; enhancer activation, including at ZGA.",
    Brd4: "BET bromodomain reader; releases paused Pol II and helps drive zygotic genome activation.",
    Chd8: "ATP-dependent chromatin remodeler.",
    Mta2: "NuRD-complex subunit (repressive chromatin remodeling).",
    Nono: "Paraspeckle RNA-binding protein; splicing and nuclear RNA retention.",
    Dpy30: "COMPASS/MLL H3K4-methylation subunit.",
    Ddx20: "DEAD-box RNA helicase (SMN/Gemin complex; snRNP assembly).",
    Rbm8a: "Exon-junction-complex core (Y14); mRNA surveillance and translation.",
    Snrpd3: "Sm-core snRNP protein; spliceosome.",
    Phf5a: "Spliceosome component (SF3b).",
    Tardbp: "TDP-43 RNA-binding protein (splicing/stability).",
    Fmr1: "FMRP translational regulator (RNA binding).",
    Msi2: "Musashi-2 RNA-binding translational regulator.",
    Tfdp1: "DP1 — E2F dimerization partner; drives cell-cycle S-phase entry.",
    Cdkn2c: "INK4 CDK inhibitor (p18); cell-cycle brake.",

    // ── sperm / paternal / germ-cell ──
    Brdt: "Testis bromodomain protein that remodels sperm chromatin — a paternal-lineage marker.",
    Ddx43: "Germ-cell DEAD-box RNA helicase (HAGE).",
    Fthl17f: "Ferritin-heavy-like X-linked gene enriched with the paternal pronucleus.",
    Nanos2: "Germ-cell RNA-binding fate determinant.",
    Hspa2: "Testis HSP70 chaperone of spermatogenesis.",
    Spz1: "Testis bHLH-Zip transcription factor.",
    Uba1y: "Y-linked ubiquitin-activating enzyme.",

    // ── other notable early-embryo genes ──
    Parp12: "PARP-family mono-ADP-ribosyltransferase (interferon/stress response).",
    Clock: "Circadian bHLH-PAS transcription factor (CLOCK/BMAL).",
    Otx1: "Anterior/head patterning homeobox transcription factor (later development).",
  };

  const desc = (g) => (g && Object.prototype.hasOwnProperty.call(DESC, g)) ? DESC[g] : null;

  // ── auto-tagger: set native title on gene-bearing elements, without clobbering existing titles ──
  function isGeneSelect(sel) {
    if (/gene/i.test(sel.id || "") || /gene/i.test(sel.name || "")) return true;
    let hit = 0; const opts = sel.options || [];
    for (let i = 0; i < opts.length && i < 40; i++) if (desc(opts[i].value)) { hit++; if (hit >= 3) return true; }
    return false;
  }
  function syncSelect(sel) {
    const apply = () => { const d = desc(sel.value); sel.title = d || ""; };
    apply();
    if (!sel._giBound) { sel._giBound = true; sel.addEventListener("change", apply); }
  }
  // tag ONE element (idempotent; never clobbers a richer existing title, e.g. a stats tooltip)
  function tagEl(el) {
    if (!el || el.nodeType !== 1 || el._giDone) return;
    if (el.hasAttribute("data-gene") || el.hasAttribute("data-g")) {
      el._giDone = 1; if (!el.title) { const d = desc(el.getAttribute("data-gene") || el.getAttribute("data-g")); if (d) el.title = d; } return;
    }
    if (el.tagName === "OPTION") { el._giDone = 1; const d = desc(el.value); if (d && !el.title) el.title = d; return; }
    if (/^(SPAN|B|STRONG|A|TD|LI)$/.test(el.tagName) && !el.children.length && !el.title) {
      el._giDone = 1; const t = el.textContent.trim(); if (t.length <= 18 && desc(t)) el.title = desc(t);
    }
  }
  function tagRoot(root) {
    if (!root || root.nodeType !== 1 && root.nodeType !== 9) return;
    if (root.nodeType === 1) tagEl(root);                                  // the node itself
    if (root.querySelectorAll) {
      root.querySelectorAll("[data-gene],[data-g],option,span,b,strong,a,td,li").forEach(tagEl);
      root.querySelectorAll("select").forEach((sel) => { if (isGeneSelect(sel)) syncSelect(sel); });
    }
    if (root.tagName === "SELECT" && isGeneSelect(root)) syncSelect(root);   // reflect selected gene on the closed control
  }

  let queued = false; const pending = [];
  function flush() { queued = false; const nodes = pending.splice(0); nodes.forEach(tagRoot); }
  function schedule(node) { pending.push(node); if (!queued) { queued = true; requestAnimationFrame(flush); } }

  function start() {
    tagRoot(document);
    const obs = new MutationObserver((muts) => {
      for (const m of muts) m.addedNodes.forEach((n) => { if (n.nodeType === 1) schedule(n); });
    });
    obs.observe(document.body, { childList: true, subtree: true });
    // options set on a select AFTER it exists (innerHTML replaced) don't always trigger for the select node;
    // a light periodic sweep of gene selects keeps their titles current cheaply.
    setInterval(() => document.querySelectorAll("select").forEach((s) => { if (isGeneSelect(s)) syncSelect(s); }), 1500);
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();

  return { desc, DESC, tag: tagRoot };
})();
