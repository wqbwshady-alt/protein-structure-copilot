/* ============================================================
   ligand.js — PDB ligand detection
   ============================================================ */

(function () {
  'use strict';

  var IGNORED_RESIDUES = new Set(['HOH', 'WAT', 'H2O', 'DOD', 'SO4', 'PO4']);

  /** Parse HETATM lines from PDB text, return sorted ligand list */
  function detectFromText(pdbText) {
    var byName = {};
    var lines = String(pdbText || '').split(/\r?\n/);
    lines.forEach(function (line) {
      if (!line || line.slice(0, 6).trim() !== 'HETATM') return;
      var resName = line.slice(17, 20).trim().toUpperCase();
      if (!resName || IGNORED_RESIDUES.has(resName)) return;
      var chainId = line.slice(21, 22).trim();
      var resId = line.slice(22, 26).trim();
      if (!byName[resName]) byName[resName] = { name: resName, count: 0, locations: [] };
      byName[resName].count += 1;
      if (chainId || resId) {
        byName[resName].locations.push((chainId ? 'chain ' + chainId + ' ' : '') + 'residue ' + resId);
      }
    });
    return Object.values(byName).sort(function (a, b) {
      return b.count - a.count || a.name.localeCompare(b.name);
    });
  }

  /** Detect ligands from a File object (returns Promise) */
  function detectFromFile(file) {
    if (!file || typeof file.text !== 'function') return Promise.resolve([]);
    return file.text().then(detectFromText).catch(function () { return []; });
  }

  /** Populate the ligand datalist + set state */
  function applyCandidates(ligands) {
    var datalist = document.getElementById('ligand-candidates');
    var ligandInput = document.getElementById('ligand-analyze');
    var names = (ligands || []).map(function (l) { return l.name; });
    var primary = names.length ? names[0] : '';

    AppState.set('ligandDetection', {
      status: primary ? 'detected' : 'none',
      ligands: ligands || [],
      primaryLigand: primary
    });

    if (datalist) {
      datalist.innerHTML = '';
      names.forEach(function (n) {
        var opt = document.createElement('option');
        opt.value = n;
        datalist.appendChild(opt);
      });
    }

    if (ligandInput && primary) ligandInput.value = primary;
    ['ligand-mutation', 'ligand-compare'].forEach(function (id) {
      var el = document.getElementById(id);
      if (el && primary) el.value = primary;
    });

    return { names: names, primary: primary };
  }

  window.LigandDetect = {
    detectFromText: detectFromText,
    detectFromFile: detectFromFile,
    applyCandidates: applyCandidates
  };
})();
