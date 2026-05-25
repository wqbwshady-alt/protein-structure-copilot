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

      h += '<div class="ranking-scroll"><div class="ranking-table">';
      h += '<div class="ranking-header">' +
        '<span>Rank</span><span>Residue</span><span>Score</span><span>Dist</span><span>Enrich</span><span>Conserv</span><span>UniProt</span>' +
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
          if (consSource === 'consurf_db') {
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

        // UniProt functional annotations
        var func = r.functional_annotations || {};
        var funcHTML = '';
        if (func.available && func.features && func.features.length > 0) {
          var featTags = func.features.slice(0, 2).map(function (f) {
            return '<span class="uniprot-tag ' + f.type.toLowerCase() + '" title="' + (f.description || '') + '">' + f.type + '</span>';
          }).join('');
          funcHTML = '<span class="rk-uniprot" title="mapping: ' + (func.mapping_confidence || 'low') + '">' + featTags + '</span>';
        } else {
          funcHTML = '<span class="rk-uniprot none" title="No UniProt functional annotation available">&mdash;</span>';
        }

        // Evidence provenance row
        var ev = r.evidence_tags || {};
        var lims = r.residue_limitations || [];
        var whyParts = [];
        if (r.why_matters) whyParts.push('<span class="ev-why">' + r.why_matters + '</span>');
        var evTags = [];
        if (ev.structural) evTags.push('<span class="ev-tag structural">[S] Structural</span>');
        if (ev.enrichment) evTags.push('<span class="ev-tag enrichment">[E] Enrichment</span>');
        if (ev.functional) evTags.push('<span class="ev-tag functional">[F] UniProt</span>');
        if (ev.conservation) evTags.push('<span class="ev-tag conservation">[C] True conservation</span>');
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
        'Key flags: static structure analysis, no energetic validation, no MD simulation, geometric classification only.' +
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

  /* ---- Init ---- */
  function init() {
    _initModeSelector();
    _initForms();
    _initDarkMode();
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
