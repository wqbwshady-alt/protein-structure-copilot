# Findings: AI Interpretation Engine Upgrade

## 1. Current Pipeline Audit

```
PDB input → parse_pdb_atoms() → ligand detection (list_ligands)
         → analyze_ligand_pocket(ligand_name, cutoff=5.0, max_interactions=15)
            → distance calculation (all ligand atoms × all protein atoms)
            → contact_residues dict (≤5.0A)
            → nearest_interactions (≤4.0A, per-residue closest contact)
            → interaction classification (classify_interaction)
            → counts by residue type
            → primary_interpretation (simple heuristic string)
         → get_hotspot_residues(interactions[:5])
         → build_report() / build_comparison_report() / build_mutation_scan_report()
         → generate_structured_interpretation(mode, data)
            → build_structured_prompt() → DeepSeek API
            → fallback: _structured_local_fallback()
         → parse_ai_sections() → OrderedDict
         → Flask route: render_index(**result) or jsonify(result)

Export: JSON (interactions array), CSV, TXT report, PyMOL script, AI report
Frontend: copilot.js renderResults() → viewer.js init() + renderHotspotList()
```

### Key Gaps Identified

| Gap | Current State | Target |
|-----|--------------|--------|
| Residue ranking | Only top-5 by distance (hotspots) | Multi-factor score per residue |
| Scoring | None — sorted by distance only | Weighted formula: distance + interaction + hotspot + diversity + chemistry |
| Confidence | "LOW/MEDIUM/HIGH" in AI text only | Structured, rule-based, per-residue + per-interaction + overall |
| Evidence citation | AI prompt asks for it, not enforced | Machine-readable evidence per interaction |
| Limitations | Generic list in AI prompt | Structured limitations dict, always present |
| AI overclaim risk | SYSTEM_PROMPT constrains, but no enforcement | Evidence-required prompt template + scoring pre-computed |

---

## 2. New Analysis Schema (v2)

### 2.1 `important_residues` Array

Each entry:
```json
{
  "chain_id": "A",
  "res_name": "ASP",
  "res_id": "25",
  "residue_key": "A:ASP25",
  "rank": 1,
  "score": 0.87,
  "score_components": {
    "distance_score": 0.28,
    "interaction_score": 0.25,
    "hotspot_score": 0.10,
    "diversity_score": 0.12,
    "chemistry_score": 0.12
  },
  "interaction_evidence": [
    {
      "type": "charged / electrostatic",
      "distance": 3.2,
      "ligand_atom": "O1",
      "protein_atom": "OD1",
      "confidence": "high"
    },
    {
      "type": "polar / possible H-bond",
      "distance": 3.5,
      "ligand_atom": "N2",
      "protein_atom": "OD2",
      "confidence": "medium"
    }
  ],
  "contact_count": 3,
  "min_distance": 3.2,
  "residue_confidence": "high",
  "pocket_location": "core",
  "why_matters": "ASP25 ranks #1 because it forms a close (3.2A) charged contact and a potential H-bond (3.5A), appearing in both hotspot and core pocket residues. Its negative charge may anchor the ligand via electrostatic steering."
}
```

### 2.2 `confidence` Block (top-level)

```json
{
  "overall_analysis_confidence": "medium",
  "confidence_reason": "PDB structure valid, ligand MK1 detected with 25 contact residues (≥10 threshold). Single static structure without MD/energy validation.",
  "residue_confidence_distribution": {
    "high": 5,
    "medium": 12,
    "low": 8
  },
  "interaction_confidence_distribution": {
    "high": 3,
    "medium": 8,
    "low": 4
  },
  "ai_interpretation_confidence": "medium",
  "data_quality_flags": {
    "has_valid_pdb": true,
    "ligand_detected": true,
    "contact_count_sufficient": true,
    "has_multiple_evidence_types": true,
    "has_energetic_validation": false,
    "has_conservation_data": false,
    "has_md_simulation": false,
    "has_multiple_structures": false
  }
}
```

### 2.3 `limitations` Block (top-level)

```json
{
  "static_structure_only": true,
  "no_energetic_validation": true,
  "no_dynamics": true,
  "no_hydrogens": true,
  "no_solvent_modeling": true,
  "no_conservation_analysis": true,
  "no_mutation_validation": true,
  "no_docking": true,
  "distance_cutoff_5A": true,
  "interaction_classification_geometric_only": true,
  "disclaimer": "This analysis is based on geometric distance criteria from a single static PDB structure. It does not constitute a validated energetic prediction. For quantitative binding assessment, MD simulation, MM-GBSA, ITC, or SPR experiments are recommended."
}
```

### 2.4 Top-level Additions (alongside existing fields)

```json
{
  "schema_version": "2.0",
  "important_residues": [...],
  "confidence": {...},
  "limitations": {...},
  "recommended_next_steps": [
    "Run 100ns MD simulation to assess pocket flexibility",
    "Perform MM-GBSA to estimate per-residue binding energy contribution",
    "Validate key contacts (ASP25, ILE50) with alanine scanning mutagenesis",
    "Compare with homologous structures to assess conservation"
  ],
  "ai_summary": "...",
  "analysis_timestamp": "2026-05-24T23:00:00Z",
  "pdb_info": {
    "filename": "1HSG.pdb",
    "source": "rcsb",
    "pdb_id": "1HSG"
  }
}
```

---

## 3. Residue Scoring Formula

### Formula

```
residue_score = distance_score + interaction_score + hotspot_score + diversity_score + chemistry_score
```

All components normalized to [0, 1]. Max possible score = 5.0.

### 3.1 distance_score (max 0.30)

```
distance_score = 0.30 × (1 - min(1, (d - d_min) / (d_max - d_min)))

Where:
  d = minimum ligand-residue contact distance (A)
  d_min = 1.5 (covalent lower bound, should use VDW radii ideally)
  d_max = 4.0 (close-contact threshold)
  Clamped to [0, 0.30]

Rationale: Closer contacts score higher. 2.0A → 0.24, 3.5A → 0.06.
```

### 3.2 interaction_score (max 0.30)

```
interaction_score = 0.30 × (weighted_contacts / max_possible_contacts)

Interaction type weights:
  charged / electrostatic  → 1.0
  polar / possible H-bond  → 0.8
  hydrophobic contact      → 0.6
  van der Waals contact    → 0.4

weighted_contacts = sum(weight_i for each contact)
max_possible_contacts = 5 × 1.0 = 5 (cap — more than 5 meaningful contacts is rare)

Rationale: Stronger interaction types contribute more. Multiple contacts increase score.
```

### 3.3 hotspot_score (max 0.15)

```
hotspot_score = 0.15 if residue is in top-5 hotspots, else 0

Rationale: Binary flag — being a hotspot means the residue already ranks high
by pure distance. This adds a small bonus to cross-validate with distance_score.
```

### 3.4 diversity_score (max 0.15)

```
diversity_score = 0.15 × unique_interaction_types / max_types

unique_interaction_types = count of distinct interaction types this residue has
max_types = 4 (charged, polar/H-bond, hydrophobic, vdW)

Rationale: A residue that participates in multiple interaction types
(e.g., both H-bond AND hydrophobic) is more important than one with
a single interaction mode. This rewards binding-site "anchor" residues.
```

### 3.5 chemistry_score (max 0.10)

```
chemistry_score = 0.10 × sum of chemistry bonuses

Chemistry bonuses (each 0.05, max 2 bonuses = 0.10):
  +0.05 if charged residue (LYS, ARG, HIS, ASP, GLU) — electrostatic steering
  +0.05 if aromatic residue (PHE, TYR, TRP) — pi-stacking potential
  +0.05 if residue is GLY in contact set — flexibility indicator
  +0.05 if residue is PRO in contact set — structural constraint

Clamped to [0, 0.10]

Rationale: Certain residue types carry intrinsic importance beyond
their interaction distance. A PHE at 3.8A with pi-stacking potential
is more significant than a GLY at 3.8A with vdW only.
```

### 3.6 Normalized Ranking

```
All residues with ≥1 interaction are scored.
Rank = descending order of residue_score.
Ties broken by: min_distance (ascending).

Score interpretation:
  0.80+  → "high importance" — likely binding anchor
  0.50–0.79 → "medium importance" — significant contributor
  0.20–0.49 → "low importance" — peripheral contact
  <0.20  → "minimal" — may be cutoff artifact
```

---

## 4. Confidence System

### 4.1 Overall Analysis Confidence

```
overall_confidence = HIGH if ALL:
  - valid PDB (is_pdb_file passes)
  - ligand detected (contact_residues not None)
  - contacts ≥ 15
  - at least 2 interaction types present
  - at least 1 charged or polar contact ≤3.5A

overall_confidence = LOW if ANY:
  - contacts < 10
  - only 1 interaction type
  - all contacts > 3.5A
  - ligand_ambiguous (multiple candidates, none dominant)

overall_confidence = MEDIUM otherwise
```

### 4.2 Per-Residue Confidence

```
residue_confidence = HIGH if ALL:
  - min_distance ≤ 3.2A
  - at least 2 contacts
  - at least 1 charged or polar/H-bond interaction

residue_confidence = LOW if ANY:
  - min_distance > 3.8A
  - only vdW contact
  - single contact only

residue_confidence = MEDIUM otherwise
```

### 4.3 Per-Interaction Confidence

```
interaction_confidence = HIGH if:
  - charged/electrostatic AND distance ≤ 3.2A
  - polar/H-bond AND distance ≤ 3.0A
  - hydrophobic AND distance ≤ 3.5A AND both atoms are carbon

interaction_confidence = LOW if:
  - vdW only AND distance > 3.5A
  - interaction classification is ambiguous (fits multiple categories)

interaction_confidence = MEDIUM otherwise
```

### 4.4 AI Interpretation Confidence

```
ai_confidence = overall_confidence, capped at MEDIUM if:
  - AI model returned from fallback (not DeepSeek)
  - AI response was truncated or malformed

ai_confidence = LOW if:
  - AI unavailable, using local rule-based fallback only
```

---

## 5. Evidence-Based AI Prompt Template

```
You are a structural biology research assistant. Your task is to interpret
pre-computed residue ranking data. You MUST cite specific evidence for
every claim.

DATA PROVIDED:
- Important residues (ranked, scored, with confidence levels)
- Interaction evidence per residue (types, distances, atoms)
- Pocket composition statistics
- Known limitations of this analysis

OUTPUT SECTIONS (use ## headers):

## Executive Summary
2-3 sentences. Name the top 3 residues with their scores and key
interaction types. State overall confidence level and primary reason.

## Top 5 Residues
For each of the top 5, write ONE sentence explaining WHY this residue
matters, citing: residue name, minimum distance, interaction type(s),
and structural role. Format: "RESNAME+ID ranks #N (score X.XX) because..."

## Evidence Assessment
Comment on the quality of structural evidence:
- Are there enough contacts to draw conclusions?
- Are the key interactions supported by multiple evidence types?
- What specific evidence is missing?

## Confidence Explanation
Explain the confidence levels. Why is the overall confidence what it is?
Which residues have high vs low confidence and why?

## Limitations
Acknowledge the specific limitations that apply to THIS analysis.
Do NOT list generic limitations — only mention ones that affect
the interpretation of THESE specific results.

## Recommended Next Steps
1-2 concrete computational or experimental follow-ups based on
the specific pocket and residues analyzed.

RULES:
1. EVERY claim about a residue MUST cite its score, min distance, and interaction type
2. If the data cannot support a conclusion, write: "Insufficient evidence to determine this"
3. NEVER invent: protein function, disease association, drug names, literature
4. Use cautious language: suggests, may, is consistent with, possibly
5. Total response: 400–500 words MAXIMUM
```

---

## 6. File-Level Implementation Plan

### Phase 1: Schema + Ranking Core

**`analysis_core.py`** — Add `ResidueRanker` class
- New class ~150 lines
- Input: `contact_residues` dict + `interactions` list
- Output: `important_residues` list with scores
- Does NOT modify existing `analyze_ligand_pocket()` or `get_hotspot_residues()`
- Preserves backward compatibility

### Phase 2: Confidence + Limitations

**`analysis_core.py`** — Add `ConfidenceAssessor` class + `LimitationsBuilder`
- New classes ~100 lines
- Input: analysis results + PDB metadata
- Output: confidence dict + limitations dict

### Phase 3: AI Prompt Integration

**`ai_client.py`** — Update `build_structured_prompt()` and `_structured_local_fallback()`
- Replace Section D "Residue-Level Evidence" with pre-computed ranking data
- Add evidence citation requirements to SYSTEM_PROMPT
- Update fallback to use new ranking data
- ~50 lines changed

### Phase 4: Report/Export Integration

**`reports.py`** — Update `build_report()`, `build_comparison_report()`, `build_mutation_scan_report()`
- Add `important_residues` section to TXT report
- Add confidence + limitations sections
- ~30 lines changed

**`app.py`** — Update `_build_analyze_result()`, `_build_protein_only_result()`, mutation_scan and compare routes
- Call new ranking/confidence/limitations functions
- Add new fields to response dicts
- ~40 lines changed

### Phase 5: Frontend Panels

**`templates/index.html`** — Add data injection for new schema fields
- `<script id="important-residues" type="application/json">...</script>`
- `<script id="confidence-data" type="application/json">...</script>`
- `<script id="limitations-data" type="application/json">...</script>`
- ~10 lines added

**`static/js/copilot.js`** — Update `_renderResults()` to show new panels
- Important Residues ranking table
- Confidence badge + summary
- Limitations collapsible card
- "Why this residue matters" on hotspot click
- ~50 lines added

**`static/css/main.css`** — New component styles
- `.confidence-badge` (HIGH=green, MEDIUM=amber, LOW=red)
- `.ranking-table` (compact monospace table)
- `.evidence-card` / `.limitations-card` styles
- ~80 lines added

### Phase 6: Tests + Code Review

**`tests/test_analysis_core.py`** — Add tests
- `test_residue_ranker_scores` — verify scoring formula
- `test_residue_ranker_ordering` — verify rank order
- `test_confidence_high` — high confidence case
- `test_confidence_low` — low confidence case (few contacts)
- `test_backward_compatibility` — old functions unchanged
- `test_ranking_with_1HSG` — real PDB integration test

**`tests/test_app.py`** — Add smoke tests
- `test_analyze_returns_new_schema_fields` — new fields in response
- `test_old_json_still_valid` — old fields still present

---

## 7. Test Plan

### Unit Tests (Phase 6)

| Test | Input | Expected |
|------|-------|----------|
| `test_distance_score` | d=2.0, d=3.5, d=5.0 | 0.24, 0.06, 0.00 |
| `test_interaction_score` | 1 charged + 1 polar | 0.30 × (1.8/5) = 0.108 |
| `test_hotspot_score` | is_hotspot=True | 0.15 |
| `test_diversity_score` | 2 unique types | 0.15 × (2/4) = 0.075 |
| `test_chemistry_score` | charged + aromatic | 0.10 |
| `test_full_score` | Mock residue with all data | Expected sum |
| `test_ranking_order` | 3 residues with known scores | Correct rank order |
| `test_confidence_high` | d=2.8, charged, 2 contacts, 20 total | "high" |
| `test_confidence_low` | d=4.5, vdW only, 1 contact, 5 total | "low" |
| `test_backward_compat` | Call analyze_ligand_pocket | Same output format as before |
| `test_old_json_fields` | Response dict keys | All old keys present |

### Integration Tests

| Test | PDB | Ligand | Check |
|------|-----|--------|-------|
| `test_1HSG_full_pipeline` | 1HSG | MK1 | Ranking + confidence + limitations present |
| `test_7VV6_full_pipeline` | 7VV6 | 6IB | Ranking + multi-ligand handling |
| `test_protein_only` | 1HSG | (none) | protein_only mode still works, ranking empty |
| `test_mutation_scan` | 1HSG | MK1, I50A | mutation fields + new schema fields |
| `test_comparison` | 1HSG WT+Mut | MK1 | comparison + new schema fields |

### Frontend Manual Checklist

- [ ] Important Residues table renders
- [ ] Confidence badge visible with correct color
- [ ] Evidence card shows per-residue interactions
- [ ] Limitations card is collapsible
- [ ] Hotspot click still focuses 3D view
- [ ] Export buttons still work (JSON includes new fields)
- [ ] 3Dmol viewer works (no regression)
- [ ] Mode switching preserves results
- [ ] Dark mode renders correctly

### Example Test Cases

**1HSG + MK1:**
- Should detect ~25 contact residues
- Top residue should score 0.70+
- At least 3 residues with confidence "high"
- Overall confidence should be "medium" or "high"
- Limitations should mention no energetic validation

**7VV6 + 6IB:**
- Should detect 6IB as ligand (CLR also present)
- Ranking should work with correct ligand
- Confidence should note multi-ligand ambiguity

**Protein-only mode:**
- `important_residues` should be empty array
- `confidence.overall` should be "low" or "insufficient_data"
- Existing protein_summary fields unchanged

---

## 8. Backward Compatibility Strategy

### Rule: New fields are additive, never remove existing fields

| Existing Field | Status |
|---------------|--------|
| `interaction_data` | Unchanged |
| `hotspot_residues` | Unchanged |
| `ai_sections` | Updated content, same structure |
| `result_text` | Unchanged |
| `report_download_url` | Unchanged |
| `json_download_url` | Updated JSON content |
| `mutation_scan_result` | Unchanged |
| `comparison_text` | Unchanged |

| New Field | Required? | Default when missing |
|-----------|-----------|---------------------|
| `important_residues` | No | `[]` |
| `confidence` | No | `null` |
| `limitations` | No | `null` |
| `recommended_next_steps` | No | `[]` |
| `ai_summary` | No | `""` |
| `schema_version` | No | `"1.0"` |

Frontend checks: `if (data.important_residues && data.important_residues.length > 0) { ... }`

---

## 9. Scientific Safety Layer (Phase 3 Detail)

### Evidence Hierarchy

```
Level 1: DIRECT STRUCTURAL EVIDENCE
  Source: PDB coordinates
  Examples: atom distances, interaction types, contact counts
  Tag: [S]
  Rule: Always present. Cannot be disputed.

Level 2: INFERRED STRUCTURAL INSIGHT
  Source: Multiple Level-1 data points combined
  Examples: pocket classification, binding mode hypothesis
  Tag: [I]
  Rule: Requires ≥2 Level-1 data points to support.

Level 3: SPECULATIVE INTERPRETATION
  Source: Level-2 insights + domain knowledge
  Examples: druggability, affinity trends, mutation impact
  Tag: [H]
  Rule: MUST include disclaimer. Cannot be stated as fact.

Level 4: FORBIDDEN CLAIMS
  Examples: numeric Kd/Ki/IC50/ΔG, disease mechanism, drug names
  Rule: AI must never make these. Post-process filter catches violations.
```

### AI Output Post-Process Validation

```
validate_ai_output(text):
  1. Forbidden phrase check:
     - /\bKd\s*[=~<>]\s*\d/ → FLAG
     - /\bIC50\b/ → FLAG
     - /\binhibits .* disease\b/ → FLAG
     - /\b(treats|cures|therapy for)\b/ → FLAG

  2. Evidence check:
     - Claims about "binding affinity" without distance citation → WARN
     - Claims about "mutation causes" without structural evidence → WARN

  3. Disclaimer check:
     - If text contains Level-3 claims and no disclaimer → APPEND disclaimer

  4. Tag application (fallback output only):
     - Prefix each paragraph with [S], [I], or [H] based on content
```

### Unsafe Output Examples (for testing)

| Unsafe Phrase | Why Unsafe | Replacement |
|--------------|-----------|-------------|
| "ASP25 binds MK1 with Kd ~ 50 nM" | Numeric affinity without measurement | "ASP25 forms close contacts (3.2A) with MK1, suggesting binding importance" |
| "This mutation causes drug resistance" | Disease claim | "This mutation alters the local physicochemical environment, which may affect ligand binding" |
| "The pocket is identical to HIV protease" | Unverified comparison | "The pocket composition (primarily hydrophobic, 25 contacts) is consistent with a small-molecule binding site" |
| "Ritonavir would bind here" | Drug name | No replacement — never name drugs |

---

## 10. Frontend State Architecture (Phase 6 Detail)

### AnalysisResult Type

```typescript
// Conceptual type — actual implementation is plain JS objects
AnalysisResult = {
  success: bool,
  analysis_mode: 'ligand' | 'protein_only' | 'mutation' | 'comparison',
  ligand_name: string,
  schema_version: '2.0',

  // Legacy fields (unchanged)
  pdb_url: string | null,
  result_text: string,
  ai_html: string,
  interaction_data: Interaction[],
  hotspot_residues: HotspotResidue[],
  lost_residues: ResidueRef[],
  gained_residues: ResidueRef[],
  ai_sections: { [title: string]: string },
  report_download_url: string | null,
  json_download_url: string | null,
  csv_download_url: string | null,
  ai_report_download_url: string | null,

  // Mutation-specific (unchanged)
  mutation_scan_result: object | null,
  mutation_scan_text: string | null,

  // Comparison-specific (unchanged)
  comparison_text: string | null,

  // Protein-only (unchanged)
  protein_summary: object | null,

  // v2 NEW FIELDS
  important_residues: ImportantResidue[],
  confidence: ConfidenceBlock | null,
  limitations: LimitationsBlock | null,
  recommended_next_steps: string[],
  ai_summary: string
}
```

### View State Transitions

```
IDLE ──[user clicks Analyze]──▶ LOADING ──[API responds]──▶ RESULT
  │                                │                          │
  │                                │ [error]                  │ [user switches mode]
  │                                ▼                          │
  │                              ERROR                       ▼
  │                                                        IDLE (new mode)
  │                                                          │
  │                                                          │ [cached result exists]
  │                                                          ▼
  └──────────────────────────────────────────────────── RESULT (restored)
```

### Viewer Lifecycle

```
create:
  PSCViewer.init({ pdbUrl, ligandName, ... })
  → viewer.clear() if existing
  → fetch(pdbUrl) → $3Dmol.createViewer('viewer', ...)
  → viewer.addModel(pdbData)
  → redraw()

update:
  Toolbar buttons call redraw() which re-applies all styles
  No re-fetch of PDB data (data is in 3Dmol memory)

destroy:
  viewer.clear() called before new init
  No explicit destroy needed (3Dmol handles cleanup)

regression prevention:
  - viewer.js only exports init(), redraw(), initFromServer()
  - No external code can modify viewer state directly
  - All style functions are private (_applyBaseStyle, _drawLines, etc.)
```

---

## 11. Compatibility Gate (Phase 8 Detail)

### Old JSON Compatibility

```python
# Test: old JSON can be parsed by new frontend
old_json = {
    "success": True,
    "interaction_data": [...],
    "hotspot_residues": [...],
    # No important_residues, confidence, limitations, schema_version
}
# New frontend: if (!data.important_residues) → skip new panels
# Old frontend: ignores unknown keys → works fine
```

### Schema Version Detection

```javascript
// copilot.js
function _renderResults(data) {
    // Always render old sections (hotspots, AI, exports)
    _renderLegacyResults(data);

    // Conditionally render new v2 panels
    if (data.schema_version === '2.0' || data.important_residues) {
        _renderImportantResiduesPanel(data.important_residues);
        _renderConfidenceBadge(data.confidence);
        _renderLimitationsCard(data.limitations);
    }
}
```

---

## 12. Governance Rules

### Code Review Gate Checklist

```
[ ] No numeric affinity values in any output path
[ ] No disease/drug/literature claims
[ ] All AI output includes disclaimer
[ ] All new JSON fields are optional (frontend handles missing)
[ ] Old JSON fields unchanged (interaction_data, hotspot_residues)
[ ] No existing function signatures modified
[ ] All new code uses Cryo-EM Dark CSS variables
[ ] No inline styles in new HTML
[ ] JS modules have single responsibility
[ ] No circular dependencies between JS modules
```

### YAGNI Enforcement

```
NOT building:
  ✗ Conservation scoring (no MSA pipeline yet)
  ✗ Energy calculation (FoldX/Rosetta not integrated)
  ✗ Multi-structure alignment (single-structure MVP)
  ✗ User preferences storage (no accounts)
  ✗ Analysis history beyond recent list (stats.json is enough)
  ✗ Real-time collaboration
  ✗ Mobile responsive beyond current @media query
  ✗ i18n / localization
  ✗ Dark/light theme beyond toggle (dark is default)
  ✗ Plugin system for scoring components
```

---

## 13. Test Data

### Test PDB: 1HSG + MK1
- 25 contact residues, 10 close interactions (≤4.0A)
- Mixed pocket: hydrophobic + polar + charged
- Known hotspot: ILE50, ILE84, ASP25, GLY27, ALA28
- Good for: scoring formula validation, confidence HIGH scenario

### Test PDB: 7VV6 + 6IB
- Multiple ligands (6IB, CLR)
- Good for: ligand ambiguity handling, multi-ligand detection

### Test PDB: protein-only (1HSG, skip ligand)
- Good for: protein_only mode regression, empty ranking

### Edge Cases
- PDB with no HETATM → list_ligands returns [], skip-ligand auto-enabled
- PDB with only water → same as above
- PDB with 100+ contacts → scoring formula still works, ranking still sorted
- PDB with 1 contact → confidence LOW, ranking has 1 entry
- PDB file > 10MB → upload warning, analysis still works
