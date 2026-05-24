# Task Plan: AI-Native Scientific Platform Upgrade

**Status:** Architecture phase — awaiting user approval
**Created:** 2026-05-24 · **Updated:** 2026-05-24
**Goal:** Upgrade "contact listing + AI vague interpretation" to "evidence-based scientific platform with scoring, confidence, safety guardrails, and modular frontend state"

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────┐
│                    FRONTEND LAYER                        │
│  ┌──────────┐  ┌──────────┐  ┌──────────────────────┐  │
│  │ Viewer   │  │ Copilot  │  │ Important Residues    │  │
│  │ (3Dmol)  │  │ (Results)│  │ Panel + Confidence    │  │
│  └────┬─────┘  └────┬─────┘  └──────────┬───────────┘  │
│       │              │                   │              │
│  ┌────┴──────────────┴───────────────────┴──────────┐   │
│  │              AppState (state.js)                  │   │
│  └────────────────────────┬─────────────────────────┘   │
├───────────────────────────┼─────────────────────────────┤
│                    API LAYER                             │
│  ┌────────────────────────┴─────────────────────────┐   │
│  │              api.js (unified fetch)               │   │
│  └────────────────────────┬─────────────────────────┘   │
├───────────────────────────┼─────────────────────────────┤
│                  FLASK BACKEND                           │
│  ┌──────────┐  ┌──────────┐  ┌──────────────────────┐  │
│  │ app.py   │  │ reports  │  │ ai_client.py          │  │
│  │ (routes) │  │ (export) │  │ (prompt + fallback)   │  │
│  └────┬─────┘  └────┬─────┘  └──────────┬───────────┘  │
│       │              │                   │              │
│  ┌────┴──────────────┴───────────────────┴──────────┐   │
│  │           analysis_core.py                        │   │
│  │  ┌──────────┐ ┌────────────┐ ┌────────────────┐  │   │
│  │  │ Residue  │ │ Confidence │ │ Safety         │  │   │
│  │  │ Ranker   │ │ Assessor   │ │ Guardrails     │  │   │
│  │  └──────────┘ └────────────┘ └────────────────┘  │   │
│  └──────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
```

---

## Phases

| # | Phase | Risk | Files | Status |
|---|-------|------|-------|--------|
| 1 | Schema + Ranking Core | MEDIUM | analysis_core.py | **complete** |
| 2 | Confidence + Limitations | LOW | analysis_core.py | **complete** |
| 3 | Scientific Safety Layer | MEDIUM | analysis_core.py, ai_client.py | **complete** |
| 4 | AI Prompt Integration | LOW | ai_client.py | **complete** |
| 5 | Report/Export Integration | LOW | reports.py, app.py | **complete** |
| 6 | Frontend State Refactor | MEDIUM | state.js, copilot.js, viewer.js, api.js, upload.js | **complete** |
| 7 | Frontend UI Panels | MEDIUM | index.html, copilot.js, main.css | **complete** |
| 8 | Compatibility + Regression | HIGH | tests/*.py, all JS files | **complete** |
| 9 | Code Review + Stabilize | LOW | all modified files | **complete** |

### Completion Definition per Phase

Each phase is **complete** when:
1. Code passes existing tests (no regression)
2. New tests for that phase pass
3. Flask server starts without errors
4. Manual smoke test: 1HSG + MK1 analysis returns correct new fields
5. Code review gate: no scientific overclaiming, no broken backward compat

### Rollback Strategy

Each phase adds **new** code only. No existing function signatures change. If a phase introduces a problem:
- Phase 1-2: comment out `ResidueRanker` / `ConfidenceAssessor` calls in app.py
- Phase 3: safety guardrails are additive, toggling off leaves original prompts intact
- Phase 4-5: new response fields are `if`-guarded in frontend
- Phase 6: state.js is standalone, reverting copilot.js to previous version restores old behavior
- Phase 7: UI panels are `if (data.important_residues)` gated
- Phase 8: tests are new files, existing tests unchanged

### Files Modified (complete)

| File | Phases | Risk Level | Rollback |
|------|--------|-----------|----------|
| analysis_core.py | 1, 2, 3 | MEDIUM | New classes only; old functions unchanged |
| ai_client.py | 3, 4 | LOW | Prompt templates updated; fallbacks preserved |
| reports.py | 5 | LOW | Add sections to report text; old sections kept |
| app.py | 5 | LOW | Add new keys to response dicts; old keys kept |
| templates/index.html | 7 | MEDIUM | Add new panels; existing panels unchanged |
| static/js/state.js | 6 | MEDIUM | Add new state keys |
| static/js/copilot.js | 6, 7 | MEDIUM | Add render methods; existing methods kept |
| static/js/viewer.js | 6 | LOW | Add event listener; existing init unchanged |
| static/js/api.js | 6 | LOW | May add helper; existing calls unchanged |
| static/css/main.css | 7 | LOW | New class additions only |
| tests/test_analysis_core.py | 8 | NONE | New test file |
| tests/test_app.py | 8 | NONE | Add test cases |
| tests/test_regression.py | 8 | NONE | New test file |

---

## A. Backward Compatibility Gate

### Phase 8: Compatibility + Regression

**Goal:** Guarantee no existing functionality silently breaks.

### Regression Checklist

```
FORM INPUTS:
  [ ] Local PDB upload (.pdb file)
  [ ] RCSB PDB ID fetch (4-char ID)
  [ ] Ligand name input (text + datalist)
  [ ] Skip ligand checkbox
  [ ] Mutation text input (e.g. R273H)
  [ ] Chain ID input
  [ ] WT + Mutant dual input

ANALYSIS MODES:
  [ ] Single structure + ligand (1HSG + MK1)
  [ ] Single structure protein-only (1HSG, skip ligand)
  [ ] Mutation scan (1HSG + MK1 + I50A)
  [ ] WT vs Mutant comparison (2 structures + ligand)

3D VIEWER:
  [ ] Protein cartoon renders (spectrum)
  [ ] Ligand renders as red sticks
  [ ] Pocket residues render as yellow sticks
  [ ] Interaction lines draw (yellow/blue/lime)
  [ ] Hotspot residues render orange
  [ ] Distance labels show
  [ ] Reset View button
  [ ] Focus Ligand button
  [ ] Focus Hotspots button
  [ ] Show/Hide Surface button
  [ ] Lost/gained residue highlights (comparison mode)

AI + RESULTS:
  [ ] AI sections render (7 sections)
  [ ] Hotspot list renders clickable
  [ ] AI fallback works when DeepSeek unavailable
  [ ] Protein-only result renders

EXPORT:
  [ ] Report TXT download
  [ ] JSON download (old fields present)
  [ ] CSV download
  [ ] AI Report download
  [ ] PyMOL script generation

API:
  [ ] /health returns 200
  [ ] /api/stats returns counts
  [ ] /api/recent_analyses returns list
  [ ] /analyze POST (server-rendered)
  [ ] /analyze_async POST (AJAX)
  [ ] /mutation_scan POST
  [ ] /compare POST
  [ ] /fetch_pdb POST

UI:
  [ ] Mode switching (3 modes)
  [ ] Drag-and-drop file upload
  [ ] File card shows after upload
  [ ] Remove file button
  [ ] Progress overlay during analysis
  [ ] Error display
  [ ] Dark mode toggle
  [ ] Recent analyses bar
  [ ] Stats counter
```

### Compatibility Strategy

1. **Schema versioning**: response includes `"schema_version": "2.0"` when new fields present. Frontend checks this before rendering new panels. Old frontend ignores unknown keys.
2. **Additive-only**: no existing response keys are removed, renamed, or changed in type.
3. **Feature flags**: new analysis functions are called via wrapper that catches exceptions and logs without breaking the main pipeline.
4. **Fallback behavior**: if `ResidueRanker` raises, `important_residues = []`. If `ConfidenceAssessor` raises, `confidence = null`. Old fields (`hotspot_residues`, `interaction_data`) are always computed first and never depend on new code.

---

## B. Scientific Safety Layer

### Phase 3: Scientific Safety + Interpretation Guardrails

**Goal:** Prevent AI from making scientifically unsupported claims.

### Safety Architecture

```
USER INPUT (PDB + ligand)
        │
        ▼
┌───────────────────────────┐
│  Deterministic Analysis    │  ← analysis_core.py (NO LLM)
│  - parse_pdb_atoms         │
│  - analyze_ligand_pocket   │
│  - ResidueRanker.score()   │
│  - ConfidenceAssessor()    │
│  - SafetyGuardrails.check()│
└───────────┬───────────────┘
            │ pre-computed data only
            ▼
┌───────────────────────────┐
│  AI Interpretation         │  ← ai_client.py (LLM with guardrails)
│  - build_safe_prompt()     │
│  - enforce_evidence_rules()│
│  - validate_ai_output()    │
│  - attach_disclaimer()     │
└───────────────────────────┘
```

### Evidence Hierarchy

```
Level 1: DIRECT STRUCTURAL EVIDENCE (always present)
  - Atom coordinates, distances, interaction types
  - Contact counts, residue composition
  → Can say: "ASP25 forms a charged contact at 3.2A"

Level 2: INFERRED STRUCTURAL INSIGHT (requires multiple data points)
  - Pocket classification (hydrophobic/charged/mixed)
  - Binding mode hypothesis (shape complementarity, electrostatic steering)
  - Flexibility assessment (Gly/Pro content)
  → Can say: "The pocket suggests hydrophobic-driven binding"

Level 3: SPECULATIVE INTERPRETATION (needs external validation disclaimer)
  - Druggability assessment
  - Affinity prediction
  - Mutation impact on binding
  - Comparison to known protein families
  → MUST say: "This is a structural hypothesis. Energetic validation is required."

Level 4: FORBIDDEN CLAIMS (AI must never make)
  - Specific binding affinity values (e.g. "Kd = 50 nM")
  - Disease mechanism assertions
  - Drug efficacy claims
  - Literature citations not in training data
  - Comparisons to unverified structures
```

### Reasoning Constraints (injected into every AI prompt)

```
CONSTRAINTS — VIOLATIONS WILL BE FLAGGED:
1. Every claim MUST cite at least one of: residue name+number, distance(A), interaction type
2. If data is insufficient for a conclusion, write: "Insufficient structural evidence to determine this"
3. NEVER state a numeric binding affinity, Kd, Ki, IC50, or ΔG value
4. NEVER claim a mutation causes/prevents disease
5. NEVER name a specific drug or therapeutic candidate
6. NEVER cite scientific literature, papers, or databases
7. ALWAYS distinguish between "structural observation" and "functional interpretation"
8. ALWAYS include the phrase "based on static structural analysis" when making inferences
9. When discussing mutation impact, ALWAYS state: "This is a heuristic assessment, not a validated ΔΔG prediction"
10. End EVERY interpretation with the limitations disclaimer
```

### Heuristic vs Validated Marking System

Every output block is tagged:

| Tag | Meaning | Example |
|-----|---------|---------|
| `[S]` | Structural observation (fact from PDB) | "ASP25 is 3.2A from ligand atom O1" |
| `[I]` | Inference (reasoning from structure) | "This suggests electrostatic steering" |
| `[H]` | Heuristic (rule-based prediction) | "Mutation likely weakens binding" |
| `[E]` | Would require experimental validation | "SPR or ITC needed to confirm affinity" |

### Hallucination Prevention

1. **Temperature = 0.3** for DeepSeek calls (lower creativity)
2. **Pre-compute everything the LLM could hallucinate** — ranking, confidence, limitations are all deterministic before the prompt is built
3. **Post-process AI output**: check for forbidden phrases (Kd, IC50, disease names, drug names) and strip/warn
4. **Fallback always available**: `_structured_local_fallback()` uses zero LLM, pure template filling
5. **AI output is labeled**: `[AI-generated interpretation — verify before citing]`

### Example: Safe vs Unsafe Output

**UNSAFE (overclaiming):**
> "ASP25 is a critical catalytic residue that binds the ligand with high affinity (estimated Kd ~ 10 nM). This is consistent with HIV protease inhibitors like Darunavir."

**SAFE (evidence-grounded):**
> "ASP25 forms a charged contact with ligand atom O1 at 3.2A and a potential H-bond with ligand atom N2 at 3.5A [S]. This suggests ASP25 may serve as an electrostatic anchor in the binding pocket [I]. Whether this translates to high binding affinity cannot be determined from structural data alone [E]. This is a structural hypothesis, not a validated energetic prediction."

---

## C. Frontend State Refactor

### Phase 6: Frontend State Refactor

**Goal:** Prevent ligand disappearing, white screens, stale state, duplicate renders.

### State Architecture

```
AppState (state.js) — SINGLE SOURCE OF TRUTH
│
├── analysis: { ready, name, source, fileObject, pdbFilename, pdbId }
├── wt: { ready, name, source, fileObject, pdbFilename, pdbId }
├── mutant: { ready, name, source, fileObject, pdbFilename, pdbId }
├── mode: 'single' | 'mutation' | 'compare'
├── ligandDetection: { status, ligands[], primaryLigand }
│
├── results: {
│     single: AnalysisResult | null,
│     mutation: AnalysisResult | null,
│     compare: AnalysisResult | null
│   }
│
├── viewer: {
│     needsInit: bool,
│     pdbUrl: string | null,
│     ligandName: string,
│     interactions: [],
│     hotspots: [],
│     lostResidues: [],
│     gainedResidues: [],
│     importantResidues: [],
│     confidence: null,
│     limitations: null
│   }
│
├── ui: {
│     loading: bool,
│     error: string | null,
│     activePanel: 'input' | 'results',
│     skipLigand: bool
│   }
│
└── stats: { totalAnalyses: 0, recent: [] }
```

### State Ownership Diagram

```
upload.js ──WRITES──▶ AppState.analysis
                      AppState.wt
                      AppState.mutant
                      AppState.ligandDetection

copilot.js ─WRITES──▶ AppState.mode
                      AppState.results[mode]
                      AppState.viewer
                      AppState.ui
                      AppState.stats

viewer.js ──READS───▶ AppState.viewer
           ─WRITES──▶ (none — viewer is render-only)

api.js ──READS──────▶ AppState.analysis (for building FormData)
        ─WRITES─────▶ (none — api.js is pure fetch)
```

### Event Flow

```
User clicks "Analyze"
  → copilot.js: AppState.set('ui', { loading: true })
  → copilot.js: _buildFormData() reads AppState.analysis/wt/mutant
  → api.js: API.analyze(formData)
  → copilot.js: receives response data
  → copilot.js: AppState.set('viewer', { pdbUrl, ligandName, interactions, ... })
  → copilot.js: AppState.set('results.' + mode, data)
  → copilot.js: AppState.set('ui', { loading: false })
  → viewer.js (via AppState.on): detects viewer.pdbUrl change
  → viewer.js: PSCViewer.init({ pdbUrl, ligandName, ... })
  → viewer.js: _renderHotspotList() reads AppState.viewer.hotspots
```

### Render Lifecycle

```
1. copilot._renderResults(data)
   ├── Builds HTML from data (AI sections, hotspots, mutation, comparison, export btns)
   ├── AppState.set('viewer', { pdbUrl, ligandName, interactions, ... })
   └── PSCViewer.init({ pdbUrl, ligandName, ... })

2. PSCViewer.init(opts)
   ├── Clears previous viewer (viewer.clear())
   ├── Fetches PDB file
   ├── $3Dmol.createViewer('viewer', ...)
   ├── viewer.addModel(pdbData, 'pdb')
   ├── redraw() → applyBaseStyle → drawLines → drawLabels → drawHotspots → drawMutations
   └── _renderHotspotList() → DOM buttons

3. AppState.on('viewer') — NOT USED for re-render (prevents loops)
   Viewer re-render only happens via explicit PSCViewer.init() call
```

### Async Strategy

```
Analysis flow:
  [User Click] → showProgress() → API.analyze() → progress.done() → renderResults() → PSCViewer.init()

PDB fetch flow:
  [User Click Fetch] → setStatus('Fetching...') → API.fetchPDB() → onSuccess(state) → setStatus('Loaded')

Stats flow:
  [Page Load] → API.fetchStats() + API.fetchRecent() → DOM update
  [After Analysis] → same refresh
  [Every 60s] → same refresh (silent)

All API calls:
  - 60s timeout (AbortController)
  - Error → user-visible message
  - No silent failures
```

### Prevented Bugs

| Bug | Prevention |
|-----|-----------|
| Ligand disappearing after mode switch | AppState.results[mode] preserves viewer state; _restoreResults re-inits viewer |
| Hotspot list empty after re-render | _renderHotspotList called AFTER DOM update in renderResults |
| Double 3Dmol init | viewer.js checks `if (viewer) { viewer.clear(); viewer = null; }` before creating new |
| Stale analysis state | AppState.set('viewer', ...) always called before PSCViewer.init() |
| White screen on error | showError() always called in .catch(); error box has display:block fallback |
| Duplicate event listeners | All addEventListener calls are in init functions called once on page load |

---

## D. Updated Phase Structure

### Phase 1: Schema + Ranking Core

**Goal:** Implement ResidueRanker with deterministic scoring.

**Risk:** MEDIUM — new logic in core analysis file.
**Files:** analysis_core.py (~150 new lines)
**Compatibility:** New class, does not modify existing functions. Called via wrapper in app.py.
**Rollback:** Remove ResidueRanker call from app.py; old hotspot logic unchanged.
**Complete when:**
- `ResidueRanker.score(contact_residues, interactions)` returns sorted list
- Each residue has: rank, score, score_components, interaction_evidence
- Unit tests pass for all 5 score components
- 1HSG + MK1 integration test returns valid ranking
**Test:** `test_residue_ranker.py` — 7 unit tests + 1 integration test

### Phase 2: Confidence + Limitations

**Goal:** Implement ConfidenceAssessor + LimitationsBuilder.

**Risk:** LOW — rule-based, no external dependencies.
**Files:** analysis_core.py (~100 new lines)
**Compatibility:** New classes. Called after ranking.
**Rollback:** Remove calls from app.py; confidence defaults to null.
**Complete when:**
- Overall, per-residue, per-interaction, AI confidence all computed
- Limitations dict built with accurate flags
- Unit tests pass for HIGH/MEDIUM/LOW scenarios
**Test:** `test_confidence.py` — 5 unit tests

### Phase 3: Scientific Safety Layer

**Goal:** Implement SafetyGuardrails + evidence hierarchy.

**Risk:** MEDIUM — changes how AI prompts are built.
**Files:** analysis_core.py (~80 lines), ai_client.py (~40 lines)
**Compatibility:** Guardrails are injected into existing prompt template. Fallback unchanged.
**Rollback:** Toggle off guardrails → original prompt behavior restored.
**Complete when:**
- `SafetyGuardrails.validate_ai_output()` catches forbidden claims
- Evidence hierarchy [S][I][H][E] tags applied to fallback output
- Prompt constraints list injected into all AI calls
- Example safe/unsafe outputs documented
**Test:** `test_safety.py` — 4 unit tests (forbidden phrase detection, tag application, disclaimer presence)

### Phase 4: AI Prompt Integration

**Goal:** Update AI prompts to consume ranking data + safety constraints.

**Risk:** LOW — template updates.
**Files:** ai_client.py (~60 lines changed)
**Compatibility:** Old prompt sections preserved. New sections use pre-computed data.
**Rollback:** Revert prompt template; old sections still functional.
**Complete when:**
- `build_structured_prompt()` includes ranking table in data block
- Section D "Residue-Level Evidence" uses pre-computed ranking instead of raw interactions
- `_structured_local_fallback()` uses ranking data
- Fallback output includes [S][I][H][E] tags
- DeepSeek prompt includes safety constraints
**Test:** Integration test — verify AI prompt contains ranking data and constraints

### Phase 5: Report/Export Integration

**Goal:** Add new fields to reports and JSON exports.

**Risk:** LOW — additive changes only.
**Files:** reports.py (~40 lines), app.py (~50 lines)
**Compatibility:** Old report sections preserved. New sections appended.
**Rollback:** New fields are optional in frontend; removing them falls back to old display.
**Complete when:**
- TXT report includes "Important Residues (Ranked)" section
- TXT report includes "Confidence Assessment" section
- TXT report includes "Limitations" section
- JSON export includes `important_residues`, `confidence`, `limitations`, `schema_version`
- Old JSON fields (interaction_data, hotspot_residues) unchanged
- All Flask routes return new fields in response
**Test:** `test_export.py` — verify JSON structure, verify old fields present

### Phase 6: Frontend State Refactor

**Goal:** Centralize state management, prevent viewer regressions.

**Risk:** MEDIUM — touches multiple JS files.
**Files:** state.js (~30 new keys), copilot.js (~60 lines changed), viewer.js (~20 lines changed), api.js (~10 lines), upload.js (~20 lines)
**Compatibility:** Existing AppState keys unchanged. New keys are additive.
**Rollback:** Revert JS files to previous versions; old behavior restored.
**Complete when:**
- AppState has all keys from section C
- Viewer init reads from AppState.viewer (not scattered globals)
- Mode switch preserves viewer state via AppState.results[mode]
- No duplicate 3Dmol init
- No white screens on error
**Test:** Manual checklist — mode switch, error state, double-analyze, page reload

### Phase 7: Frontend UI Panels

**Goal:** Render Important Residues panel, confidence badges, limitations card.

**Risk:** MEDIUM — new DOM elements, new CSS.
**Files:** templates/index.html (~15 lines), copilot.js (~80 lines), main.css (~80 lines)
**Compatibility:** New panels inside `#results-area`, gated by `if (data.important_residues)`.
**Rollback:** Remove panel HTML; old results display unaffected.
**Complete when:**
- Important Residues ranking table renders
- Confidence badge with color (HIGH=#00e699, MEDIUM=#f0a020, LOW=#f0475b)
- Evidence card shows per-residue interaction details
- Limitations card is collapsible
- "Why this residue matters" appears on hotspot hover/click
- All panels work in dark mode
- No 3Dmol viewer regression
**Test:** Manual visual checklist — 9 items

### Phase 8: Compatibility + Regression Tests

**Goal:** Prove nothing is broken.

**Risk:** HIGH — this is the gate before claiming "done".
**Files:** tests/test_regression.py (new), tests/test_analysis_core.py (expanded), tests/test_app.py (expanded)
**Compatibility:** New tests only; existing tests run first to confirm no regression.
**Rollback:** N/A — this phase only adds tests.
**Complete when:**
- All 30+ regression checklist items verified
- Old JSON format still parses correctly
- All 3 analysis modes produce valid output
- All export formats include both old and new fields
- Frontend renders correctly with and without new fields
**Test:** `test_regression.py` — 30+ automated checks

### Phase 9: Code Review + Stabilization

**Goal:** Final quality gate.

**Risk:** LOW — review only, no new features.
**Files:** All modified files reviewed.
**Compatibility:** Review checklist includes backward compat verification.
**Rollback:** Fix issues found, re-run regression tests.
**Complete when:**
- Code review passes (no HIGH severity issues)
- No scientific overclaiming in any output path
- No duplicate code across JS modules
- All tests pass (existing + new)
- Schema v2 documented
- Manual smoke test passes

---

## E. Engineering Governance

### Code Review Gates

| Gate | When | What |
|------|------|------|
| Pre-commit | Each phase | Run existing tests, verify no regression |
| Phase review | End of phase | Review diff, verify phase completion definition |
| Safety review | Phases 3, 4 | Check AI output for overclaiming, forbidden phrases |
| Compatibility review | Phase 8 | Full regression checklist |
| Final review | Phase 9 | Complete code review + stabilization |

### No-Overengineering Rules

1. Three similar lines > one premature abstraction
2. No interface/abstract class until 3+ implementations exist
3. No config file until 5+ configurable values
4. No plugin system until 3+ plugins exist
5. No microservice split until single server hits real limits
6. Scoring formula: 5 components, not 15 — resist the urge to add "one more factor"

### YAGNI Checks (per phase)

Before implementing, ask:
- Does the current user need this? (not a hypothetical future user)
- Can it be added later without breaking anything?
- Is there a simpler version that works for now?

### Schema Freeze Policy

- Schema v2 is FROZEN after Phase 5 completion
- New fields added as optional (v2.1, v2.2) — never remove or rename
- Breaking changes require v3 with migration path
- Schema version in every JSON response

### UI Consistency Rules

- All new components use Cryo-EM Dark variables (--accent, --data, --bg-elevated, etc.)
- No inline styles (use CSS classes)
- All text uses DM Sans (UI) or JetBrains Mono (data)
- Confidence colors: HIGH=--data (#00e699), MEDIUM=--warning (#f0a020), LOW=--error (#f0475b)
- All panels collapsible if > 5 items

### AI Output Quality Rules

- Every claim cites specific structural evidence
- No numeric affinity values
- No disease/drug/literature claims
- Disclaimer always present
- Fallback output clearly labeled as non-AI

### Testing Requirements

- Unit tests: 1 per scoring component, 1 per confidence level
- Integration tests: 1 per analysis mode (single, mutation, comparison, protein-only)
- Regression tests: full checklist (30+ items)
- Manual smoke test: 1HSG + MK1 before marking any phase complete

### Rollback Requirements

- Each phase must be independently revertible
- New code is additive (new classes, new functions) — not modifying existing signatures
- Feature flag: `AppState.get('features').rankingEnabled` to toggle new panels

---

## F. Final Deliverables

### 1. Updated Architecture Roadmap

```
NOW (Phase 1-3):         Core engine (ranking + confidence + safety)
    ↓
NEXT (Phase 4-5):        Integration (AI prompts + exports)
    ↓
THEN (Phase 6-7):        Frontend (state refactor + UI panels)
    ↓
GATE (Phase 8-9):        Quality (regression + code review)
    ↓
FUTURE (Post-v2):        Multi-structure alignment, conservation scoring,
                         FoldX/Rosetta integration, MD trajectory viewer,
                         Database of analyzed structures, Public deployment
```

### 2. Updated Implementation Phases

See Section D above — 9 phases, each with risk level, files, compatibility strategy, rollback, completion definition, and test strategy.

### 3. Technical Debt Map

| Debt | Severity | Location | Fix Phase |
|------|----------|----------|-----------|
| Duplicate `resolve_loaded_structure` / `resolve_structure_input` | LOW | app.py:123-159 | Future |
| stats.js still exists but not loaded | LOW | static/js/stats.js | Delete in Phase 6 |
| AI response can contain JSON-breaking control chars | MEDIUM | ai_client.py | Future |
| No request rate limiting | MEDIUM | app.py | Future |
| results/ directory grows unbounded | LOW | app.py | Future |
| No PDB file cleanup (uploads/ grows) | LOW | app.py | Future |
| index.html inline data scripts are verbose | LOW | index.html | Phase 7 |

### 4. Highest-Risk Areas

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Scoring formula breaks old hotspot ordering | Low | High | Old hotspot logic unchanged; new ranking is separate |
| AI prompt changes cause worse output | Medium | Medium | Fallback always works; prompt changes are additive |
| Frontend state refactor introduces viewer bugs | Medium | High | Phase 6 before Phase 7; test after each JS file change |
| Schema v2 fields confuse old export consumers | Low | Low | All new fields optional; schema_version signals v2 |
| Scientific overclaiming in AI output | Medium | High | Phase 3 safety guardrails; post-process AI output |

### 5. Recommended Next Immediate Step

**Phase 1: Schema + Ranking Core**

This is the foundation. Everything else depends on having scored, ranked residues. Start with `ResidueRanker` class in `analysis_core.py`, write tests, verify with 1HSG + MK1. This phase has zero frontend impact and can be merged independently.

### 6. What Should NOT Be Built Yet

- Conservation scoring (requires MSA pipeline)
- FoldX/Rosetta energy integration (requires external tools)
- MD trajectory viewer (requires MD engine)
- Multi-structure alignment viewer
- Public database of analyzed structures
- User accounts / saved analyses
- Real-time collaborative analysis
- Mobile app

### 7. MVP Boundary

```
MVP (this upgrade cycle):
  ✓ Residue ranking with deterministic scoring
  ✓ Structured confidence system
  ✓ Scientific safety guardrails
  ✓ Evidence-grounded AI interpretation
  ✓ Important Residues frontend panel
  ✓ Confidence badges + limitations card
  ✓ Backward-compatible JSON schema v2
  ✓ Regression test suite (30+ checks)

NOT in MVP:
  ✗ Conservation analysis
  ✗ Energy calculation (FoldX, MM-GBSA)
  ✗ MD simulation integration
  ✗ Multi-structure alignment
  ✗ Public database
```

### 8. Long-Term Platform Vision

```
Protein Structure Copilot v1 (current):
  "Upload PDB → see contacts → get AI interpretation"

Protein Structure Copilot v2 (this upgrade):
  "Upload PDB → ranked residues with evidence → confidence-rated AI interpretation"

Protein Structure Copilot v3 (next):
  "Upload PDB → ranked residues → confidence → conservation → mutant prediction → MD-ready"

Protein Structure Copilot v4 (vision):
  "AI-native structural biology workstation: natural language → automated pipeline →
   multi-method validation → publication-ready figures → collaboration → database"
```
