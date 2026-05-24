/* ============================================================
   upload.js — Drop zones, file selection, structure state sync
   ============================================================ */

(function () {
  'use strict';

  var MAX_FILE_SIZE = 10 * 1024 * 1024; // 10 MB warning

  function formatSize(bytes) {
    if (bytes < 1024) return bytes + ' B';
    if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + ' KB';
    return (bytes / (1024 * 1024)).toFixed(1) + ' MB';
  }

  /* ---- Apply file state to drop zone DOM ---- */
  function _applyZoneDOM(zone, structure) {
    var nameEl = zone.querySelector('.file-card-filename');
    var sizeEl = zone.querySelector('.file-card-size');
    var warningEl = zone.closest('.input-section') ? zone.parentElement.querySelector('.file-size-warning') : null;
    var errorEl = zone.parentElement.querySelector('.file-error-msg');

    if (structure.ready) {
      zone.classList.add('has-file');
      zone.classList.remove('has-error');
      if (nameEl) nameEl.textContent = structure.name || '';
      if (sizeEl) sizeEl.textContent = structure.size || '';
      if (warningEl) warningEl.style.display = 'none';
      if (errorEl) errorEl.style.display = 'none';
    } else {
      zone.classList.remove('has-file');
    }
  }

  /* ---- Update structure-source indicator ---- */
  function _updateSourceIndicator(elId, structure, label) {
    var el = document.getElementById(elId);
    if (!el) return;
    if (structure.ready) {
      el.textContent = (label || '') + ': ' + (structure.name || '') + ' loaded';
      el.className = 'status-msg success';
    } else {
      el.textContent = 'Load ' + (label || 'structure') + ' (RCSB or upload)';
      el.className = 'status-msg warning';
    }
  }

  function _setStatus(msg, type) {
    var el = document.getElementById('analysis-status');
    if (!el) return;
    el.textContent = msg;
    el.className = 'status-msg ' + (type || '');
  }

  /* ---- Build structure state from a File ---- */
  function _fileToState(file) {
    return {
      name: file.name,
      size: formatSize(file.size),
      source: 'local',
      fileObject: file,
      pdbFilename: null,
      pdbText: null,
      pdbId: null,
      ready: true
    };
  }

  /* ---- Build structure state from RCSB fetch result ---- */
  function _rcsbToState(data) {
    return {
      name: data.pdb_id + '.pdb',
      size: 'from RCSB',
      source: 'rcsb',
      fileObject: null,
      pdbFilename: data.filename,
      pdbId: data.pdb_id,
      pdbText: null,
      ready: true
    };
  }

  /* ---- Set up a single drop zone ---- */
  function _setupDropZone(zone, stateKey, onFileLoaded) {
    var input = zone.querySelector('input[type="file"]');
    var removeBtn = zone.querySelector('.remove-file-btn');

    function _handleFile(file) {
      if (!file) return;
      if (!file.name.toLowerCase().endsWith('.pdb')) {
        zone.classList.add('has-error');
        zone.classList.remove('has-file');
        var errEl = zone.parentElement.querySelector('.file-error-msg');
        if (errEl) errEl.style.display = 'block';
        return;
      }
      zone.classList.remove('has-error');
      var errEl = zone.parentElement.querySelector('.file-error-msg');
      if (errEl) errEl.style.display = 'none';

      var state = _fileToState(file);
      if (file.size > MAX_FILE_SIZE) {
        var warnEl = zone.parentElement.querySelector('.file-size-warning');
        if (warnEl) warnEl.style.display = 'block';
      }
      AppState.replace(stateKey, state);
      _applyZoneDOM(zone, state);
      if (onFileLoaded) onFileLoaded(state, file);
    }

    input.addEventListener('change', function () {
      if (input.files && input.files[0]) _handleFile(input.files[0]);
    });

    ['dragenter', 'dragover'].forEach(function (eName) {
      zone.addEventListener(eName, function (e) { e.preventDefault(); zone.classList.add('drag-over'); });
    });
    zone.addEventListener('dragleave', function (e) { e.preventDefault(); zone.classList.remove('drag-over'); });
    zone.addEventListener('drop', function (e) {
      e.preventDefault();
      zone.classList.remove('drag-over');
      if (e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0]) {
        _handleFile(e.dataTransfer.files[0]);
      }
    });

    if (removeBtn) {
      removeBtn.addEventListener('click', function (e) {
        e.stopPropagation();
        input.value = '';
        zone.classList.remove('has-file', 'has-error');
        AppState.replace(stateKey, AppState.emptyStructure());
        _applyZoneDOM(zone, AppState.emptyStructure());
        if (onFileLoaded) onFileLoaded(AppState.emptyStructure(), null);
      });
    }
  }

  /* ---- RCSB fetch handler ---- */
  function _setupRCSBFetch(btnId, inputId, statusId, stateKey, onLoaded) {
    var btn = document.getElementById(btnId);
    var input = document.getElementById(inputId);
    var statusEl = document.getElementById(statusId);
    if (!btn || !input) return;

    btn.addEventListener('click', function () {
      var pdbId = input.value.trim().toUpperCase();
      if (!pdbId || pdbId.length !== 4) {
        if (statusEl) { statusEl.textContent = 'Enter a 4-character PDB ID'; statusEl.className = 'status-msg error'; }
        return;
      }
      btn.disabled = true;
      if (statusEl) { statusEl.textContent = 'Fetching...'; statusEl.className = 'status-msg info'; }

      API.fetchPDB(pdbId).then(function (data) {
        btn.disabled = false;
        if (data.success) {
          if (statusEl) { statusEl.textContent = 'Loaded: ' + pdbId; statusEl.className = 'status-msg success'; }
          var state = _rcsbToState(data);
          AppState.replace(stateKey, state);
          var zoneId = stateKey === 'analysis' ? 'drop-zone-analyze'
            : stateKey === 'wt' ? 'drop-zone-wt' : 'drop-zone-mut';
          var zone = document.getElementById(zoneId);
          if (zone) _applyZoneDOM(zone, state);

          if (onLoaded) onLoaded(data, state);
        } else {
          if (statusEl) { statusEl.textContent = data.error_text; statusEl.className = 'status-msg error'; }
        }
      }).catch(function (err) {
        btn.disabled = false;
        if (statusEl) { statusEl.textContent = 'Network error: ' + (err.message || 'timeout'); statusEl.className = 'status-msg error'; }
      });
    });
  }

  /* ---- Skip ligand mode ---- */
  function _initSkipLigand() {
    var cb = document.getElementById('skip-ligand-mode');
    if (!cb) return;
    cb.addEventListener('change', function () {
      AppState.set('ui', { skipLigand: cb.checked });
      var ligandInput = document.getElementById('ligand-analyze');
      if (!ligandInput) return;
      if (cb.checked) {
        ligandInput.value = '';
        ligandInput.disabled = true;
        ligandInput.placeholder = 'Ligand skipped (protein-only)';
        _setStatus('Protein-only mode. No ligand required.', 'info');
      } else {
        ligandInput.disabled = false;
        ligandInput.placeholder = 'e.g. MK1 / CLR';
      }
    });
  }

  /* ---- Main init ---- */
  function init() {
    // Drop zones
    _setupDropZone(document.getElementById('drop-zone-analyze'), 'analysis', function (state, file) {
      if (state.ready && file) {
        LigandDetect.detectFromFile(file).then(function (ligands) {
          LigandDetect.applyCandidates(ligands);
          _setStatus(state.name + ' ready for analysis', 'success');
        });
      }
      _updateSourceIndicator('mutation-structure-source', state, 'Structure');
    });

    _setupDropZone(document.getElementById('drop-zone-wt'), 'wt', function (state) {
      _updateSourceIndicator('wt-structure-source', state, 'WT');
      _updateCompareReady();
    });

    _setupDropZone(document.getElementById('drop-zone-mut'), 'mutant', function (state) {
      _updateSourceIndicator('mut-structure-source', state, 'Mutant');
      _updateCompareReady();
    });

    // RCSB fetches
    _setupRCSBFetch('fetch-pdb-btn', 'pdb-id-input', 'fetch-status', 'analysis', function (data) {
      var ligands = (data.ligands || []).map(function (n) { return { name: n, count: 1 }; });
      LigandDetect.applyCandidates(ligands);
      _setStatus(data.pdb_id + ' loaded from RCSB', 'success');
      _updateSourceIndicator('mutation-structure-source', AppState.get('analysis'), 'Structure');
    });

    _setupRCSBFetch('mutation-fetch-pdb-btn', 'mutation-pdb-id-input', 'mutation-fetch-status', 'analysis', function (data) {
      _setStatus(data.pdb_id + ' loaded from RCSB', 'success');
      _updateSourceIndicator('mutation-structure-source', AppState.get('analysis'), 'Structure');
    });

    _setupRCSBFetch('wt-fetch-pdb-btn', 'wt-pdb-id-input', 'wt-fetch-status', 'wt', function () {
      _updateSourceIndicator('wt-structure-source', AppState.get('wt'), 'WT');
      _updateCompareReady();
    });

    _setupRCSBFetch('mut-fetch-pdb-btn', 'mut-pdb-id-input', 'mut-fetch-status', 'mutant', function () {
      _updateSourceIndicator('mut-structure-source', AppState.get('mutant'), 'Mutant');
      _updateCompareReady();
    });

    // Skip ligand
    _initSkipLigand();

    // React to analysis state changes
    AppState.on('analysis', function (state) {
      if (state.ready) _setStatus(state.name + ' ready for analysis', 'success');
      else _setStatus('', '');
      _updateSourceIndicator('mutation-structure-source', state, 'Structure');
    });

    AppState.on('wt', function (state) {
      _updateSourceIndicator('wt-structure-source', state, 'WT');
      _updateCompareReady();
    });

    AppState.on('mutant', function (state) {
      _updateSourceIndicator('mut-structure-source', state, 'Mutant');
      _updateCompareReady();
    });
  }

  function _updateCompareReady() {
    var wt = AppState.get('wt');
    var mut = AppState.get('mutant');
    var el = document.getElementById('compare-status');
    if (!el) return;
    if (wt.ready && mut.ready) {
      el.textContent = 'WT and mutant structures ready for comparison';
      el.className = 'status-msg success';
    } else {
      el.textContent = '';
      el.className = 'status-msg';
    }
  }

  /* ---- Export ---- */
  window.PSCInitUpload = init;
})();
