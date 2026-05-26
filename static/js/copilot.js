/* ============================================================
   copilot.js — Mode switching, form handling, progress, results
   ============================================================ */

(function () {
  'use strict';

  /* ---- Mode switching with result caching ---- */
  function _initModeSelector() {
    var options = document.querySelectorAll('.mode-option');
    var panels = {
      single: document.getElementById('panel-single'),
      mutation: document.getElementById('panel-mutation'),
      compare: document.getElementById('panel-compare')
    };

    options.forEach(function (opt) {
      opt.addEventListener('click', function () {
        var mode = opt.getAttribute('data-mode');
        var current = AppState.get('mode');
        if (mode === current) return;

        // Save current mode's results before switching
        _cacheResults(current);

        AppState.set('mode', mode);
        options.forEach(function (o) { o.classList.remove('active'); });
        opt.classList.add('active');

        Object.keys(panels).forEach(function (k) {
          if (panels[k]) panels[k].style.display = 'none';
        });
        if (panels[mode]) panels[mode].style.display = 'block';

        // Restore cached results for new mode
        _restoreResults(mode);
      });
    });
  }

  var _resultCache = { single: null, mutation: null, compare: null };

  function _cacheResults(mode) {
    var area = document.getElementById('results-area');
    if (area && area.innerHTML) {
      _resultCache[mode] = {
        html: area.innerHTML,
        viewerState: {
          ligandName: AppState.get('viewer').ligandName,
          interactions: AppState.get('viewer').interactions,
          hotspots: AppState.get('viewer').hotspots,
          lostResidues: AppState.get('viewer').lostResidues,
          gainedResidues: AppState.get('viewer').gainedResidues,
          pdbUrl: AppState.get('viewer').pdbUrl
        }
      };
    }
  }

  function _restoreResults(mode) {
    var cached = _resultCache[mode];
    var area = document.getElementById('results-area');
    var placeholder = document.getElementById('results-placeholder');
    var errorBox = document.getElementById('analysis-error');

    if (cached && cached.html) {
      if (area) area.innerHTML = cached.html;
      if (placeholder) placeholder.style.display = 'none';
      if (errorBox) errorBox.style.display = 'none';

      // Restore viewer
      var vs = cached.viewerState;
      if (vs && vs.pdbUrl && window.PSCViewer) {
        window.PSCViewer.init({
          pdbUrl: vs.pdbUrl,
          ligandName: vs.ligandName,
          interactions: vs.interactions,
          hotspots: vs.hotspots,
          lostResidues: vs.lostResidues,
          gainedResidues: vs.gainedResidues
        });
      }
    } else {
      if (area) area.innerHTML = '';
      if (placeholder) placeholder.style.display = 'block';
      if (errorBox) errorBox.style.display = 'none';
    }
  }

  /* ---- Progress overlay ---- */
  var _progressStages = [
    'Reading PDB structure', 'Detecting ligands', 'Calculating pocket contacts',
    'Generating interaction data', 'Preparing 3D visualization', 'AI structural interpretation'
  ];

  function _showProgress(title) {
    var overlay = document.getElementById('progress-overlay');
    if (overlay) overlay.classList.add('active');
    var t = document.getElementById('progress-title');
    if (t) t.textContent = title || 'Analyzing Structure';
    var stages = document.querySelectorAll('.progress-stage');
    stages.forEach(function (s) { s.className = 'progress-stage pending'; });
    var bar = document.getElementById('progress-bar-fill');
    if (bar) bar.style.width = '0%';
    var done = document.getElementById('progress-done-message');
    var err = document.getElementById('progress-error-message');
    if (done) done.classList.remove('show');
    if (err) err.classList.remove('show');
  }

  function _hideProgress() {
    var overlay = document.getElementById('progress-overlay');
    if (overlay) overlay.classList.remove('active');
  }

  function _startProgressSim() {
    _showProgress('Analyzing Structure');
    var timers = [];
    [400, 800, 1200, 1600, 2000, 2500].forEach(function (d, i) {
      timers.push(setTimeout(function () {
        var stages = document.querySelectorAll('.progress-stage');
        stages.forEach(function (s, j) {
          if (j < i) s.className = 'progress-stage done';
          else if (j === i) s.className = 'progress-stage active';
          else s.className = 'progress-stage pending';
        });
        var bar = document.getElementById('progress-bar-fill');
        if (bar) bar.style.width = Math.round((i / (_progressStages.length - 1)) * 100) + '%';
      }, d));
    });
    return {
      done: function () {
        timers.forEach(clearTimeout);
        var stages = document.querySelectorAll('.progress-stage');
        stages.forEach(function (s) { s.className = 'progress-stage done'; });
        var bar = document.getElementById('progress-bar-fill');
        if (bar) bar.style.width = '100%';
        var done = document.getElementById('progress-done-message');
        if (done) done.classList.add('show');
        setTimeout(_hideProgress, 800);
      },
      error: function (msg) {
        timers.forEach(clearTimeout);
        var err = document.getElementById('progress-error-message');
        if (err) { err.textContent = msg; err.classList.add('show'); }
        setTimeout(_hideProgress, 2000);
      }
    };
  }

  /* ---- Build FormData from mode ---- */
  function _buildFormData(form) {
    var fd = new FormData(form);
    var mode = AppState.get('mode');

    if (mode === 'single' || mode === 'mutation') {
      var af = AppState.get('analysis');
      if (af.ready) {
        if (af.source === 'rcsb') fd.set('pdb_filename', af.pdbFilename || '');
        else if (af.fileObject) fd.set('pdb_file', af.fileObject);
      }
    }

    if (mode === 'single') {
      var skip = document.getElementById('skip-ligand-mode');
      if (skip && skip.checked) fd.set('skipLigand', 'true');
    }

    if (mode === 'compare') {
      var wt = AppState.get('wt');
      var mut = AppState.get('mutant');
      if (wt.ready) {
        if (wt.source === 'rcsb') fd.set('wt_pdb_filename', wt.pdbFilename || '');
        else if (wt.fileObject) fd.set('wt_file', wt.fileObject);
      }
      if (mut.ready) {
        if (mut.source === 'rcsb') fd.set('mut_pdb_filename', mut.pdbFilename || '');
        else if (mut.fileObject) fd.set('mut_file', mut.fileObject);
      }
    }

    return fd;
  }

  /* ---- Form submission ---- */
  function _initForms() {
    var endpoints = {
      'form-single':   { fn: API.analyze,      mode: 'single' },
      'form-mutation': { fn: API.mutationScan,  mode: 'mutation' },
      'form-compare':  { fn: API.compare,       mode: 'compare' }
    };

    Object.keys(endpoints).forEach(function (formId) {
      var form = document.getElementById(formId);
      if (!form) return;
      var cfg = endpoints[formId];

      form.addEventListener('submit', function (e) {
        e.preventDefault();
        var fd = _buildFormData(form);
        var progress = _startProgressSim();

        cfg.fn(fd).then(function (data) {
          progress.done();
          if (data.success) {
            _renderResults(data, cfg.mode);
          } else {
            _showError(data.error_text || 'Analysis failed.');
          }
        }).catch(function (err) {
          progress.error('Network error');
          _showError('Request failed: ' + (err.message || 'timeout'));
        });
      });
    });
  }

  /* ---- Render results ---- */
  function _renderResults(data, mode) {
    var area = document.getElementById('results-area');
    if (!area) return;
    var errorBox = document.getElementById('analysis-error');
    if (errorBox) errorBox.style.display = 'none';
    var placeholder = document.getElementById('results-placeholder');
    if (placeholder) placeholder.style.display = 'none';

    var h = '';

    // AI sections
    if (data.ai_sections && Object.keys(data.ai_sections).length > 0) {
      Object.keys(data.ai_sections).forEach(function (key) {
        var body = data.ai_sections[key] || '';
        h += '<details class="ai-section" open><summary>' + key + '</summary>' +
          '<div class="section-body">' + body + '</div></details>';
      });
    } else if (data.ai_html) {
      h += '<div class="section-body" style="padding:0;">' + data.ai_html + '</div>';
    }

    // Protein summary (protein-only mode)
    if (data.protein_summary) {
      var ps = data.protein_summary;
      h += '<div class="mutation-summary" style="margin-top:12px;">' +
        '<p><strong>Structure Overview:</strong> ' + (ps.chain_count || 0) + ' chain(s), ' +
        (ps.residue_count || 0) + ' residues, ' + (ps.atom_count || 0) + ' atoms</p>' +
        '</div>';
    }

    // Hotspots
    if (data.hotspot_residues && data.hotspot_residues.length > 0) {
      h += '<h2 style="margin-top:16px;">Top Hotspots</h2><ul class="hotspot-list" id="hotspot-list"></ul>';
    }

    // Mutation scan
    if (data.mutation_scan_result) {
      var mr = data.mutation_scan_result;
      h += '<div class="mutation-summary" style="margin-top:12px;">' +
        '<p><strong>Mutation:</strong> ' + (mr.mutation || '') + ' chain ' + (mr.chain_id || '') + '</p>' +
        '<p>' + ((mr.original_residue || {}).res_name || '') + ((mr.original_residue || {}).res_id || '') +
        ' → ' + ((mr.mutant_residue || {}).res_name || '') + ((mr.mutant_residue || {}).res_id || '') + '</p>' +
        '</div>';
    }

    // Comparison
    if (data.comparison_text) {
      h += '<details class="ai-section" style="margin-top:12px;"><summary>Full Comparison Report</summary>' +
        '<div class="comparison-box">' + data.comparison_text.replace(/\n/g, '<br>') + '</div></details>';
    }

    // --- v2: Important Residues Ranking Table ---
    if (data.important_residues && data.important_residues.length > 0) {
      h += '<h2 style="margin-top:16px;">Important Residues</h2>';

      // Enrichment summary bar (pocket-level, shown once)
      var enrich = (data.important_residues[0] || {}).enrichment || {};
      var overallEnrich = enrich.overall_enrichment || {};
      if (overallEnrich.significant_types && overallEnrich.significant_types.length > 0) {
        h += '<div class="enrichment-bar">' +
          '<span class="enrichment-label">Pocket Enrichment (vs whole protein):</span> ' +
          overallEnrich.significant_types.map(function (s) {
            return '<span class="enrichment-tag">' + s + '</span>';
          }).join(' ') +
          '<span class="enrichment-test"> (Fisher exact, p &lt; 0.05)</span>' +
        '</div>';
      }

      // Overall annotation status bar
      var annoSummary = data.annotation_summary || {};
      if (annoSummary.status) {
        var annoBarClass = 'annotation-bar';
        var annoStatusText = '';
        var annoDetail = '';
        var statusMap = {
          success: 'UniProt annotations available',
          partial: 'UniProt annotations partially available',
          unavailable: 'UniProt: no annotations at mapped positions',
          no_mapping: 'UniProt: no PDB-to-UniProt mapping',
          failed: 'UniProt: fetch failed',
        };
        var statusClassMap = {
          success: 'anno-success',
          partial: 'anno-partial',
          unavailable: 'anno-unavailable',
          no_mapping: 'anno-unavailable',
          failed: 'anno-failed',
        };
        annoStatusText = statusMap[annoSummary.status] || 'UniProt: unknown status';
        if (annoSummary.residues_with_data !== undefined && annoSummary.total_residues) {
          var accList = (annoSummary.accessions_queried || []).join(', ') || 'none';
          annoDetail = annoSummary.residues_with_data + '/' + annoSummary.total_residues +
            ' residues annotated' + (accList ? ' (' + accList + ')' : '');
          if (annoSummary.dbref_found === false) {
            annoDetail += ' - PDB file lacks DBREF records';
          }
          if ((annoSummary.accessions_failed || []).length > 0) {
            annoDetail += ' - fetch failed for: ' + annoSummary.accessions_failed.join(', ');
          }
        }
        var annoClass = statusClassMap[annoSummary.status] || 'anno-unavailable';
        h += '<div class="annotation-bar ' + annoClass + '">' +
          '<span class="annotation-status-label">' + annoStatusText + '</span>' +
          (annoDetail ? ' <span class="annotation-detail">' + annoDetail + '</span>' : '') +
        '</div>';
      }

      // Flexibility summary card
      var flexSummary = data.flexibility || {};
      if (flexSummary.classification && flexSummary.classification !== 'no_data') {
        var flexCardClass = 'result-card';
        var ratioStr = flexSummary.flexibility_ratio ? ' (ratio: ' + flexSummary.flexibility_ratio.toFixed(2) + 'x)' : '';
        h += '<div class="result-card" style="margin-bottom:8px;">' +
          '<div class="result-card-header">' +
            '<span class="card-icon">&#x1FCCF;</span> Pocket Flexibility' +
          '</div>' +
          '<div style="font-size:12px;color:var(--text-secondary);line-height:1.6;">' +
            '<span style="font-weight:600;color:var(--text-primary);">' + (flexSummary.label || '') + '</span>' +
            ' &mdash; Pocket mean B-factor: <strong>' + (flexSummary.mean_b || '?') + '</strong>' +
            ' vs protein mean: <strong>' + (flexSummary.global_mean_b || '?') + '</strong>' + ratioStr +
          '</div>' +
          '<div style="margin-top:6px;display:flex;gap:12px;font-size:11px;flex-wrap:wrap;">';
        if (flexSummary.rigid_count > 0) {
          h += '<span style="color:#3b82f6;">&#x25CF; ' + flexSummary.rigid_count + ' rigid</span>';
        }
        if (flexSummary.flexible_count > 0) {
          h += '<span style="color:#f59e0b;">&#x25CF; ' + flexSummary.flexible_count + ' flexible</span>';
        }
        if (flexSummary.highly_flexible_count > 0) {
          h += '<span style="color:#ef4444;">&#x25CF; ' + flexSummary.highly_flexible_count + ' highly flexible</span>';
        }
        h += '</div></div>';
      }

      // Pi-stacking card
      var piData = data.pi_stacking || {};
      var piPi = piData.pi_pi_interactions || [];
      var catPi = piData.cation_pi_interactions || [];
      if (piPi.length > 0 || catPi.length > 0) {
        h += '<div class="result-card" style="margin-bottom:8px;">' +
          '<div class="result-card-header">' +
            '<span class="card-icon">&#x2B55;</span> Aromatic Interactions' +
          '</div>';
        piPi.forEach(function (p) {
          var typeLabel = {pi_pi_face_to_face:'Face-to-face', pi_pi_edge_to_face:'Edge-to-face', pi_pi_t_shaped:'T-shaped'}[p.type] || p.type;
          h += '<div style="font-size:11px;color:var(--text-secondary);margin-bottom:3px;">' +
            '&#x3C0;-&#x3C0; ' + typeLabel + ': <strong style="color:var(--text-primary);">' +
            p.residue1 + ' &harr; ' + p.residue2 + '</strong> ' +
            '(' + p.distance + 'A, ' + p.angle + '&deg;)' +
          '</div>';
        });
        catPi.forEach(function (c) {
          h += '<div style="font-size:11px;color:var(--text-secondary);margin-bottom:3px;">' +
            'Cation-&#x3C0;: <strong style="color:var(--text-primary);">' +
            c.cationic_residue + ' &rarr; ' + c.aromatic_residue + '</strong> ' +
            '(' + c.distance + 'A)' +
          '</div>';
        });
        h += '</div>';
      }

      // Ligand profile card
      var lp = data.ligand_profile || {};
      if (lp.mw) {
        h += '<div class="result-card" style="margin-bottom:8px;">' +
          '<div class="result-card-header">' +
            '<span class="card-icon">&#x1F9EA;</span> Ligand Profile (' + (lp.name || '') + ')' +
          '</div>' +
          '<div style="font-size:11px;color:var(--text-secondary);line-height:1.8;display:grid;grid-template-columns:1fr 1fr;gap:2px 16px;">' +
            '<span>MW: <strong>' + (lp.mw || '?') + ' Da</strong></span>' +
            '<span>LogP: <strong>' + (lp.logp || '?') + '</strong></span>' +
            '<span>TPSA: <strong>' + (lp.tpsa || '?') + ' A²</strong></span>' +
            '<span>HBD: <strong>' + (lp.hbd || 0) + '</strong> / HBA: <strong>' + (lp.hba || 0) + '</strong></span>' +
            '<span>Rotatable bonds: <strong>' + (lp.rotatable_bonds || 0) + '</strong></span>' +
            '<span>Rings: <strong>' + (lp.ring_count || 0) + '</strong> (arom: ' + (lp.aromatic_rings || 0) + ')</span>' +
            '<span>Drug-likeness: <strong>' + (lp.drug_likeness || '?') + '</strong></span>' +
            '<span>Ro5: <strong>' + ((lp.ro5_violations === 0) ? '0 violations (pass)' : (lp.ro5_violations || '?') + ' violations') + '</strong></span>' +
          '</div>';
        if (lp.mmff_strain_energy !== undefined && lp.mmff_strain_energy !== null) {
          h += '<div style="margin-top:6px;font-size:11px;color:var(--text-secondary);">' +
            'MMFF94 strain energy: <strong>' + (lp.mmff_strain_energy > 0 ? '+' : '') + lp.mmff_strain_energy.toFixed(1) + ' kcal/mol</strong>' +
            ' (bound vs relaxed conformer)' +
          '</div>';
        }
        h += '</div>';
      }

      // SASA card
      var sasa = data.sasa || {};
      if (sasa.total_sasa > 0) {
        h += '<div class="result-card" style="margin-bottom:8px;">' +
          '<div class="result-card-header">' +
            '<span class="card-icon">&#x1F4D0;</span> Solvent Accessibility (SASA)' +
          '</div>' +
          '<div style="font-size:12px;color:var(--text-secondary);line-height:1.6;">' +
            'Pocket total SASA: <strong>' + (sasa.total_sasa || 0).toFixed(0) + ' A²</strong>' +
          '</div>' +
          '<div style="margin-top:4px;font-size:10px;color:var(--text-muted);">' +
            'Shrake-Rupley algorithm, 1.4A water probe. Values relative to Gly-X-Gly tripeptide reference.' +
          '</div></div>';
      }

      // Water bridge card
      var wb = data.water_bridges || {};
      if (wb.total > 0) {
        h += '<div class="result-card" style="margin-bottom:8px;">' +
          '<div class="result-card-header">' +
            '<span class="card-icon">&#x1F4A7;</span> Water-Mediated Contacts' +
          '</div>' +
          '<div style="font-size:12px;color:var(--text-secondary);line-height:1.6;">' +
            'Total bridges: <strong>' + wb.total + '</strong>' +
            ' &mdash; Water molecules: <strong>' + (wb.water_molecules_involved || 0) + '</strong>' +
            ' &middot; Residues bridged: <strong>' + (wb.protein_residues_bridged || 0) + '</strong>' +
          '</div>' +
          '<div style="margin-top:4px;font-size:10px;color:var(--text-muted);">' +
            'Water O within 3.5A of both protein and ligand atoms. Water-mediated H-bonds can contribute significantly to binding affinity.' +
          '</div></div>';
      }

      // Salt bridge card
      var sb = data.salt_bridges || {};
      if (sb.total > 0) {
        h += '<div class="result-card" style="margin-bottom:8px;">' +
          '<div class="result-card-header">' +
            '<span class="card-icon">&#x1F9F2;</span> Salt Bridges (ASP/GLU ↔ LYS/ARG/HIS)' +
          '</div>' +
          '<div style="font-size:12px;color:var(--text-secondary);line-height:1.6;">' +
            'Total: <strong>' + sb.total + '</strong>' +
            ' &mdash; Strong: <strong style="color:#3b82f6;">' + (sb.strong || 0) + '</strong> (<3.2A)' +
            ' &middot; Moderate: <strong style="color:#8b5cf6;">' + (sb.moderate || 0) + '</strong> (3.2-3.6A)' +
            ' &middot; Weak: <strong style="color:#94a3b8;">' + (sb.weak || 0) + '</strong> (3.6-4.0A)' +
          '</div>' +
          '<div style="margin-top:4px;font-size:10px;color:var(--text-muted);">' +
            'Close-range charge-charge interactions between carboxylate (ASP/GLU) and ammonium/guanidinium (LYS/ARG/HIS) side chains.' +
          '</div></div>';
      }

      // Hydrogen bond geometry card
      var hb = data.hbonds || {};
      if (hb.total > 0) {
        h += '<div class="result-card" style="margin-bottom:8px;">' +
          '<div class="result-card-header">' +
            '<span class="card-icon">&#x1F4CC;</span> Hydrogen Bonds (Baker-Hubbard)' +
          '</div>' +
          '<div style="font-size:12px;color:var(--text-secondary);line-height:1.6;">' +
            'Total: <strong>' + hb.total + '</strong>' +
            ' &mdash; Validated: <strong style="color:#10b981;">' + (hb.validated || 0) + '</strong>' +
            ' &middot; Possible: <strong style="color:#f59e0b;">' + (hb.possible || 0) + '</strong>' +
            ' &middot; Protein-Ligand: <strong>' + (hb.protein_ligand || 0) + '</strong>' +
          '</div>' +
          '<div style="margin-top:4px;font-size:10px;color:var(--text-muted);">' +
            'H···A &lt; 2.5A, D···A &lt; 3.5A, D-H···A &gt; 120° (validated) / &gt; 90° (possible). H positions estimated from heavy-atom geometry.' +
          '</div></div>';
      }

      // Prodigy affinity card
      var prod = data.prodigy || {};
      if (prod.delta_g !== undefined || prod.kd !== undefined) {
        h += '<div class="result-card" style="margin-bottom:8px;">' +
          '<div class="result-card-header">' +
            '<span class="card-icon">&#x1F3AF;</span> Binding Affinity Prediction (Prodigy)' +
          '</div>' +
          '<div style="font-size:12px;color:var(--text-secondary);line-height:1.6;">';
        if (prod.delta_g !== undefined && prod.delta_g !== null) {
          h += 'ΔG = <strong style="color:var(--text-primary);">' + prod.delta_g + ' kcal/mol</strong>';
        }
        if (prod.kd !== undefined && prod.kd !== null) {
          h += (prod.delta_g !== undefined ? ' &mdash; ' : '') + 'Kd ≈ <strong style="color:var(--text-primary);">' + prod.kd + '</strong>';
        }
        h += '</div>' +
          '<div style="margin-top:4px;font-size:10px;color:var(--text-muted);">' +
            'Predicted by Prodigy (Utrecht University). Statistical model — not experimental. Use as qualitative reference.' +
          '</div></div>';
      }

      // Energy summary card
      var eSummary = data.interaction_energy || {};
      if (eSummary.total_energy !== undefined) {
        h += '<div class="result-card" style="margin-bottom:8px;">' +
          '<div class="result-card-header">' +
            '<span class="card-icon">&#x26A1;</span> Interaction Energy (LJ + Coulomb)' +
          '</div>' +
          '<div style="font-size:12px;color:var(--text-secondary);line-height:1.6;">' +
            'Total: <strong style="color:var(--text-primary);">' + (eSummary.total_energy || 0).toFixed(1) + ' kcal/mol</strong>' +
            ' &mdash; vdW: ' + (eSummary.total_vdw || 0).toFixed(1) +
            ', Coulomb: ' + (eSummary.total_coulomb || 0).toFixed(1) +
          '</div>' +
          '<div style="margin-top:4px;font-size:10px;color:var(--text-muted);">' +
            'Approximate gas-phase energy. Simplied AMBER-like parameters, distance-dependent dielectric (ε=4r). Use for relative ranking.' +
          '</div>' +
        '</div>';
      }

      h += '<div class="ranking-scroll"><div class="ranking-table">';
      h += '<div class="ranking-header">' +
        '<span>Rank</span><span>Residue</span><span>Score</span><span>Dist</span><span>Energy</span><span>Flex</span><span>Enrich</span><span>Conserv</span><span>UniProt</span>' +
      '</div>';
      data.important_residues.slice(0, 10).forEach(function (r) {
        var confClass = 'conf-' + (r.residue_confidence || 'low');
        var types = [];
        (r.interaction_evidence || []).forEach(function (e) { if (types.indexOf(e.type) === -1) types.push(e.type); });

        // Enrichment
        var catEnrich = (r.enrichment || {}).category_enrichment || null;
        var enrichHTML = '';
        if (catEnrich) {
          var fold = catEnrich.fold_enrichment;
          var sig = catEnrich.significant;
          var pv = catEnrich.p_value;
          if (fold !== undefined && fold > 0) {
            var enrichClass = sig ? 'enrich-sig' : 'enrich-ns';
            var enrichTitle = 'Fold: ' + fold.toFixed(1) + 'x, p=' + (pv !== null ? pv : 'N/A') +
              ', ' + (r.enrichment.category || '') + ' category';
            enrichHTML = '<span class="rk-enrich ' + enrichClass + '" title="' + enrichTitle + '">' +
              (fold >= 1 ? '+' : '') + fold.toFixed(1) + 'x</span>';
          }
        }

        // Energy
        var energy = r.interaction_energy || {};
        var energyHTML = '';
        if (energy.total !== undefined && energy.total !== 0) {
          var eClass = energy.total < -2 ? 'energy-strong' : energy.total < -0.5 ? 'energy-mod' : 'energy-weak';
          var eTitle = 'vdW: ' + (energy.vdw||0) + ' + Coulomb: ' + (energy.coulomb||0) + ' kcal/mol';
          energyHTML = '<span class="rk-energy ' + eClass + '" title="' + eTitle + '">' +
            (energy.total > 0 ? '+' : '') + energy.total.toFixed(1) + '</span>';
        } else {
          energyHTML = '<span class="rk-energy energy-none">—</span>';
        }

        // Flexibility
        var flex = r.flexibility || {};
        var flexHTML = '';
        if (flex.classification && flex.classification !== 'unknown') {
          var flexClass = 'flex-' + flex.classification;
          var flexLabel = {'rigid':'[R]','normal':'','flexible':'[F]','highly_flexible':'[FF]'}[flex.classification] || '';
          flexHTML = '<span class="rk-flex ' + flexClass + '" title="B-factor: ' + (flex.mean_b || '?') + ' (z=' + (flex.z_score || '?') + ')">' +
            (flex.mean_b ? flex.mean_b.toFixed(0) : '?') + ' <span class="flex-src">' + flexLabel + '</span></span>';
        }

        // Conservation
        var cons = r.conservation || {};
        var consHTML = '';
        if (cons.score !== undefined) {
          var consSource = cons.source || 'blosum62_proxy';
          var consAvailable = cons.available;
          var consTitle = cons.source_detail || '';
          var consScore = cons.score;
          var consClass = consAvailable ? 'cons-real' : 'cons-proxy';
          var consLabel = '';
          if (consSource === 'uniprot_variant_constraint') {
            consLabel = '[V]';
            consClass = 'cons-real cons-variant';
          } else if (consSource === 'consurf_db') {
            consLabel = '[C]';
            consClass = 'cons-real cons-consurf';
          } else if (consSource === 'blosum62_proxy') {
            consLabel = '[P]';
            consClass = 'cons-proxy';
          } else {
            consLabel = '(' + consSource + ')';
          }
          consHTML = '<span class="rk-cons ' + consClass + '" title="' + consTitle + '">' +
            consScore.toFixed(2) + ' <span class="cons-src">' + consLabel + '</span></span>';
        }

        // UniProt functional annotations — status-aware rendering
        var func = r.functional_annotations || {};
        var annoStatus = r.annotation_status || 'failed';
        var funcHTML = '';
        if (func.available && func.features && func.features.length > 0) {
          var featTags = func.features.slice(0, 2).map(function (f) {
            return '<span class="uniprot-tag ' + f.type.toLowerCase() + '" title="' + (f.description || '') + '">' + f.type + '</span>';
          }).join('');
          var statusTitle = 'mapping: ' + (func.mapping_confidence || 'low');
          funcHTML = '<span class="rk-uniprot" title="' + statusTitle + '">' + featTags + '</span>';
        } else {
          var statusLabel = {no_mapping:'no mapping', unavailable:'no annotation', failed:'fetch failed', skipped:'skipped'}[annoStatus] || 'no data';
          var statusTitle = {no_mapping:'No PDB-to-UniProt mapping', unavailable:'No annotation at mapped position', failed:'UniProt API fetch failed', skipped:'Non-standard residue skipped'}[annoStatus] || 'No UniProt data';
          funcHTML = '<span class="rk-uniprot none status-' + annoStatus + '" title="' + statusTitle + '">&mdash; <span class="uniprot-status-text">' + statusLabel + '</span></span>';
        }

        // Evidence provenance row
        var ev = r.evidence_tags || {};
        var lims = r.residue_limitations || [];
        var whyParts = [];
        if (r.why_matters) whyParts.push('<span class="ev-why">' + r.why_matters + '</span>');
        var evTags = [];
        if (ev.structural) evTags.push('<span class="ev-tag structural">[S] Structural</span>');
        if (ev.enrichment) evTags.push('<span class="ev-tag enrichment">[E] Enrichment</span>');
        if (ev.functional) {
          var funcTagClass = 'functional';
          if (annoStatus === 'partial') funcTagClass = 'functional-partial';
          evTags.push('<span class="ev-tag ' + funcTagClass + '">[F] UniProt</span>');
        } else if (annoStatus === 'failed') {
          evTags.push('<span class="ev-tag functional-failed">[F] UniProt failed</span>');
        }
        if (ev.conservation) {
          var consSrc = cons.source || '';
          var consTagLabel = consSrc === 'uniprot_variant_constraint' ? '[V] Functional constraint' : '[C] True conservation';
          evTags.push('<span class="ev-tag conservation">' + consTagLabel + '</span>');
        }
        if (ev.proxy_only) evTags.push('<span class="ev-tag proxy">[P] Substitution proxy</span>');
        if (evTags.length) whyParts.push('<span class="ev-tags">' + evTags.join(' ') + '</span>');
        if (lims.length > 0) {
          whyParts.push('<span class="ev-lims">' + lims.slice(0, 2).join('; ') + '</span>');
        }

        h += '<div class="ranking-row" data-chain="' + r.chain_id + '" data-resi="' + r.res_id + '" data-resname="' + r.res_name + '">' +
          '<span class="rk-rank">#' + r.rank + '</span>' +
          '<span class="rk-residue">' + (r.residue_key || '') + '</span>' +
          '<span class="rk-score">' + (r.score || 0).toFixed(2) + '</span>' +
          '<span class="rk-dist">' + (r.min_distance || '?') + 'A</span>' +
          energyHTML +
          flexHTML +
          enrichHTML +
          consHTML +
          funcHTML +
          (whyParts.length ? '<div class="ranking-why">' + whyParts.join(' &middot; ') + '</div>' : '') +
        '</div>';
      });
      h += '</div></div>';
    }

    // --- v2: Confidence Badge ---
    if (data.confidence) {
      var overall = data.confidence.overall_analysis_confidence || 'unknown';
      var confColor = overall === 'high' ? 'var(--data)' : overall === 'medium' ? 'var(--warning)' : 'var(--error)';
      h += '<div class="confidence-badge" style="margin-top:12px;">' +
        '<span class="conf-dot" style="background:' + confColor + ';"></span>' +
        '<span>Analysis Confidence: <strong>' + overall.toUpperCase() + '</strong></span>' +
        '<span class="conf-reason">' + (data.confidence.confidence_reason || '') + '</span>' +
      '</div>';
    }

    // --- v2: Limitations Card ---
    if (data.limitations) {
      h += '<details class="ai-section" style="margin-top:8px;">' +
        '<summary>Limitations & Caveats</summary>' +
        '<div class="section-body">' +
        '<p>' + (data.limitations.disclaimer || '') + '</p>' +
        '<p style="margin-top:8px;font-size:11px;color:var(--text-muted);">' +
        'Key flags: single static PDB structure, simplified LJ+Coulomb energy scoring (qualitative ranking), no MD simulation, no experimental binding validation.' +
        '</p></div></details>';
    }

    // Export buttons
    h += '<div class="export-group">';
    if (data.report_download_url) h += '<a href="' + data.report_download_url + '" class="export-btn export-txt">Report (TXT)</a>';
    if (data.json_download_url) h += '<a href="' + data.json_download_url + '" class="export-btn export-json">Data (JSON)</a>';
    if (data.csv_download_url) h += '<a href="' + data.csv_download_url + '" class="export-btn export-csv">Data (CSV)</a>';
    if (data.ai_report_download_url) h += '<a href="' + data.ai_report_download_url + '" class="export-btn export-pml">AI Report</a>';
    h += '</div>';

    area.innerHTML = h;
    area.classList.add('animate-in');

    // Update AppState viewer data
    AppState.set('viewer', {
      ligandName: data.ligand_name || '',
      interactions: data.interaction_data || [],
      hotspots: data.hotspot_residues || [],
      lostResidues: data.lost_residues || [],
      gainedResidues: data.gained_residues || [],
      pdbUrl: data.pdb_url || null
    });

    // Init viewer
    if (data.pdb_url && window.PSCViewer) {
      window.PSCViewer.init({
        pdbUrl: data.pdb_url,
        ligandName: data.ligand_name || '',
        interactions: data.interaction_data || [],
        hotspots: data.hotspot_residues || [],
        lostResidues: data.lost_residues || [],
        gainedResidues: data.gained_residues || []
      });
    }

    // Refresh stats
    _refreshStats();
  }

  function _showError(msg) {
    var box = document.getElementById('analysis-error');
    if (box) { box.textContent = msg; box.style.display = 'block'; }
    AppState.set('ui', { error: msg });
  }

  /* ---- Stats refresh ---- */
  function _refreshStats() {
    API.fetchStats().then(function (data) {
      var el = document.getElementById('stats-count');
      if (el) _animateCount(el, data.total_analyses || 0);
    }).catch(function () {});

    API.fetchRecent().then(function (data) {
      var list = document.getElementById('recent-list');
      if (!list) return;
      var items = (data || []).slice(0, 8);
      if (!items.length) { list.innerHTML = '<span style="color:var(--text-muted);font-size:12px;">No analyses yet</span>'; return; }
      list.innerHTML = items.map(function (item) {
        var tc = { single: 'single', mutation: 'mutation', comparison: 'comparison' };
        var tl = { single: 'ANALYSIS', mutation: 'MUTATION', comparison: 'COMPARE' };
        var t = item.analysis_type || 'single';
        var pid = (item.pdb_id || item.pdb_name || '');
        if (pid.length > 24) pid = pid.slice(0, 22) + '..';
        return '<div class="recent-item">' +
          '<span class="recent-type ' + (tc[t] || 'single') + '">' + (tl[t] || 'ANALYSIS') + '</span>' +
          '<span class="recent-pdb">' + pid + '</span>' +
          '<span class="recent-time">' + _timeAgo(item.timestamp) + '</span></div>';
      }).join('');
    }).catch(function () {});
  }

  function _animateCount(el, target) {
    var cur = parseInt(el.textContent, 10) || 0;
    if (cur === target) { el.textContent = target; return; }
    var step = Math.max(1, Math.floor(Math.abs(target - cur) / 20));
    var iv = setInterval(function () {
      if (cur < target) cur = Math.min(cur + step, target);
      else cur = Math.max(cur - step, target);
      el.textContent = cur;
      if (cur === target) clearInterval(iv);
    }, 30);
  }

  function _timeAgo(iso) {
    if (!iso) return '';
    var diff = Math.floor((Date.now() - new Date(iso).getTime()) / 1000);
    if (diff < 60) return 'just now';
    if (diff < 3600) return Math.floor(diff / 60) + 'm ago';
    if (diff < 86400) return Math.floor(diff / 3600) + 'h ago';
    return Math.floor(diff / 86400) + 'd ago';
  }

  /* ---- Dark mode ---- */
  function _initDarkMode() {
    var toggle = document.getElementById('dark-mode-toggle');
    if (toggle) toggle.addEventListener('click', function () { document.body.classList.toggle('dark'); });
  }

  /* ---- Demo button ---- */
  function _initDemoButton() {
    var btn = document.getElementById('demo-btn');
    if (!btn) return;
    btn.addEventListener('click', function () {
      // Fill PDB ID and ligand
      var pdbInput = document.getElementById('pdb-id-input');
      var ligandInput = document.getElementById('ligand-analyze');
      if (pdbInput) pdbInput.value = '1ATP';
      if (ligandInput) ligandInput.value = 'ATP';

      btn.textContent = 'Loading example...';
      btn.disabled = true;

      // Fetch 1ATP from RCSB
      var fetchBtn = document.getElementById('fetch-pdb-btn');
      if (fetchBtn) {
        fetchBtn.click();
      }

      // Wait for fetch to complete, then analyze
      var checkReady = setInterval(function () {
        var af = AppState.get('analysis');
        if (af && af.ready) {
          clearInterval(checkReady);
          var form = document.getElementById('form-single');
          if (form) {
            var submitEvent = new Event('submit', { bubbles: true, cancelable: true });
            form.dispatchEvent(submitEvent);
          }
          btn.textContent = 'Try with Example';
          btn.disabled = false;
        }
      }, 300);

      // Timeout after 15s
      setTimeout(function () {
        clearInterval(checkReady);
        btn.textContent = 'Try with Example';
        btn.disabled = false;
      }, 15000);
    });
  }

  /* ---- Init ---- */
  function init() {
    _initModeSelector();
    _initForms();
    _initDarkMode();
    _initDemoButton();
    _refreshStats();

    // Server-rendered results
    if (window.PSCViewer && window.PSCViewer.initFromServer) {
      window.PSCViewer.initFromServer();
    }

    // Periodic stats refresh
    setInterval(_refreshStats, 60000);
  }

  window.PSCInitCopilot = init;
})();
