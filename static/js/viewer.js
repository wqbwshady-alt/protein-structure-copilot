/* ============================================================
   viewer.js — 3Dmol.js viewer management

   Public API:
     PSCViewer.init(opts)  — { pdbUrl, ligandName, interactions, hotspots, lostResidues, gainedResidues }
     PSCViewer.initFromServer() — reads <script id="pdb-url-data"> etc. from DOM
     PSCViewer.redraw()
     PSCViewer.focusResidue(residue)
   ============================================================ */

(function () {
  'use strict';

  var viewer = null;
  var surfaceOn = false;
  var _pocketSelections = new Set();

  // Current render state
  var _ligandName = '';
  var _interactions = [];
  var _hotspots = [];
  var _lostResidues = [];
  var _gainedResidues = [];

  /* ---- Helpers ---- */
  function _ixColor(ix) { return ix.color || 'yellow'; }

  /* ---- Base styling ---- */
  function _applyBaseStyle() {
    if (!viewer) return;
    viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
    if (!_ligandName) return;

    // Ligand: bright red sticks
    viewer.setStyle({ resn: _ligandName }, { stick: { color: 'red', radius: 0.25 } });

    // Pocket residues: yellow sticks within 5A
    viewer.addStyle(
      { within: { distance: 5, sel: { resn: _ligandName } } },
      { stick: { color: 'yellow', radius: 0.18 } }
    );

    _pocketSelections.clear();
    _interactions.forEach(function (ix) {
      _pocketSelections.add(ix.chain_id + ':' + ix.res_id);
    });
  }

  /* ---- Draw interaction lines ---- */
  function _drawLines() {
    if (!viewer) return;
    _interactions.forEach(function (ix) {
      viewer.addCylinder({
        start: { x: ix.start[0], y: ix.start[1], z: ix.start[2] },
        end:   { x: ix.end[0],   y: ix.end[1],   z: ix.end[2] },
        radius: 0.08, color: _ixColor(ix), dashed: true
      });
    });
  }

  /* ---- Draw distance labels ---- */
  function _drawLabels() {
    if (!viewer) return;
    var top5 = _hotspots.slice(0, 5);
    var shown = 0;
    _interactions.forEach(function (ix) {
      if (shown >= 5) return;
      if (!top5.some(function (r) { return r.chain_id === ix.chain_id && String(r.res_id) === String(ix.res_id); })) return;
      shown++;
      var mx = (ix.start[0] + ix.end[0]) / 2;
      var my = (ix.start[1] + ix.end[1]) / 2;
      var mz = (ix.start[2] + ix.end[2]) / 2;
      viewer.addLabel(ix.distance.toFixed(1) + ' A', {
        position: { x: mx, y: my, z: mz }, fontColor: _ixColor(ix), fontSize: 10, showBackground: false
      });
    });
  }

  /* ---- Draw hotspots ---- */
  function _drawHotspots() {
    if (!viewer) return;
    _hotspots.forEach(function (r) {
      viewer.setStyle(
        { chain: r.chain_id, resi: r.res_id },
        { stick: { color: 'orange', radius: 0.45 }, sphere: { color: 'orange', radius: 0.75 } }
      );
      viewer.addResLabels(
        { chain: r.chain_id, resi: r.res_id },
        { fontColor: 'black', backgroundColor: 'orange', fontSize: 12, showBackground: true }
      );
    });
  }

  /* ---- Draw mutation highlights ---- */
  function _drawMutations() {
    if (!viewer) return;
    _lostResidues.forEach(function (r) {
      viewer.setStyle({ chain: r.chain_id, resi: r.res_id }, { stick: { color: 'red', radius: 0.35 } });
    });
    _gainedResidues.forEach(function (r) {
      viewer.setStyle({ chain: r.chain_id, resi: r.res_id }, { stick: { color: 'green', radius: 0.35 } });
    });
  }

  /* ---- Pocket surface ---- */
  function _showSurface() {
    if (!viewer) return;
    _pocketSelections.forEach(function (key) {
      var parts = key.split(':');
      viewer.addSurface($3Dmol.SurfaceType.VDW, { opacity: 0.4, color: 'yellow' }, { chain: parts[0], resi: parts[1] });
    });
  }

  /* ---- Hotspot list ---- */
  function _renderHotspotList() {
    var list = document.getElementById('hotspot-list');
    if (!list) return;
    list.innerHTML = '';
    if (!_hotspots || _hotspots.length === 0) {
      list.innerHTML = '<li style="color:var(--text-muted);font-size:13px;padding:8px;">No hotspot residues detected.</li>';
      return;
    }
    _hotspots.forEach(function (r, i) {
      var li = document.createElement('li');
      var btn = document.createElement('button');
      btn.className = 'hotspot-btn';
      btn.innerHTML = '<span class="hotspot-rank">#' + (i + 1) + '</span>' +
        r.res_name + r.res_id + ' chain ' + r.chain_id + ' · ' +
        r.distance + ' A · ' + (r.interaction_type || 'contact');
      btn.addEventListener('click', function () {
        _focusResidue(r);
      });
      li.appendChild(btn);
      list.appendChild(li);
    });
  }

  /* ---- Full redraw ---- */
  function redraw() {
    if (!viewer) return;
    viewer.removeAllShapes();
    viewer.removeAllLabels();
    viewer.removeAllSurfaces();

    _applyBaseStyle();
    _drawLines();
    _drawLabels();
    _drawHotspots();
    _drawMutations();

    if (surfaceOn) _showSurface();

    if (_ligandName) viewer.zoomTo({ resn: _ligandName });
    else viewer.zoomTo();
    viewer.render();
    viewer.resize();

    var sb = document.getElementById('surface-toggle-btn');
    if (sb) {
      sb.textContent = surfaceOn ? 'Hide Surface' : 'Show Surface';
      sb.classList.toggle('active', surfaceOn);
    }
  }

  /* ---- Focus single residue ---- */
  function _focusResidue(residue) {
    if (!viewer) return;
    viewer.removeAllShapes();
    viewer.removeAllLabels();
    viewer.removeAllSurfaces();

    viewer.setStyle({}, { cartoon: { color: 'lightgrey', opacity: 0.32 } });
    viewer.setStyle({ resn: _ligandName }, { stick: { color: 'yellow', radius: 0.3 } });
    viewer.setStyle(
      { chain: residue.chain_id, resi: residue.res_id },
      { stick: { color: 'orange', radius: 0.55 }, sphere: { color: 'orange', radius: 1.0 } }
    );
    viewer.addResLabels(
      { chain: residue.chain_id, resi: residue.res_id },
      { fontColor: 'black', backgroundColor: 'orange', fontSize: 14, showBackground: true }
    );

    _interactions.filter(function (ix) {
      return ix.chain_id === residue.chain_id && String(ix.res_id) === String(residue.res_id);
    }).forEach(function (ix) {
      viewer.addCylinder({
        start: { x: ix.start[0], y: ix.start[1], z: ix.start[2] },
        end:   { x: ix.end[0],   y: ix.end[1],   z: ix.end[2] },
        radius: 0.13, color: _ixColor(ix), dashed: true
      });
      var mx = (ix.start[0] + ix.end[0]) / 2;
      var my = (ix.start[1] + ix.end[1]) / 2;
      var mz = (ix.start[2] + ix.end[2]) / 2;
      viewer.addLabel(ix.distance.toFixed(1) + ' A', {
        position: { x: mx, y: my, z: mz }, fontColor: _ixColor(ix), fontSize: 11, showBackground: false
      });
    });

    viewer.zoomTo({ chain: residue.chain_id, resi: residue.res_id });
    viewer.render();
  }

  /* ---- Focus all hotspots ---- */
  function _focusHotspotsView() {
    if (!viewer) return;
    viewer.removeAllShapes();
    viewer.removeAllLabels();
    viewer.removeAllSurfaces();

    viewer.setStyle({}, { cartoon: { color: 'lightgrey', opacity: 0.25 } });
    viewer.setStyle({ resn: _ligandName }, { stick: { color: 'yellow', radius: 0.3 } });
    _drawHotspots();
    _drawLines();
    _drawLabels();

    if (_hotspots.length > 0) {
      viewer.zoomTo(_hotspots.map(function (r) { return { chain: r.chain_id, resi: r.res_id }; }));
    } else if (_ligandName) {
      viewer.zoomTo({ resn: _ligandName });
    }
    viewer.render();
  }

  /* ---- Init viewer with PDB data ---- */
  function init(opts) {
    opts = opts || {};
    _ligandName = opts.ligandName || '';
    _interactions = opts.interactions || [];
    _hotspots = opts.hotspots || [];
    _lostResidues = opts.lostResidues || [];
    _gainedResidues = opts.gainedResidues || [];

    if (!opts.pdbUrl) return;

    if (viewer) { viewer.clear(); viewer = null; }

    fetch(opts.pdbUrl)
      .then(function (r) { return r.text(); })
      .then(function (pdbData) {
        viewer = $3Dmol.createViewer('viewer', { backgroundColor: '#080d1a' });
        viewer.addModel(pdbData, 'pdb');

        var emptyEl = document.getElementById('viewer-empty-state');
        if (emptyEl) emptyEl.style.display = 'none';

        redraw();
        _renderHotspotList();
      })
      .catch(function (err) {
        console.error('Failed to load PDB:', err);
      });
  }

  /* ---- Init from server-rendered data (replaces index.html inline script) ---- */
  function initFromServer() {
    var pdbUrlEl = document.getElementById('pdb-url-data');
    if (!pdbUrlEl) return;
    var pdbUrl = pdbUrlEl.getAttribute('data-url');
    if (!pdbUrl) return;

    var ligandName = (document.getElementById('ligand-name-data') || {}).getAttribute('data-value') || '';
    var interactions = JSON.parse((document.getElementById('interaction-data') || { textContent: '[]' }).textContent);
    var lostResidues = JSON.parse((document.getElementById('lost-residues') || { textContent: '[]' }).textContent);
    var gainedResidues = JSON.parse((document.getElementById('gained-residues') || { textContent: '[]' }).textContent);
    var hotspots = JSON.parse((document.getElementById('hotspot-residues') || { textContent: '[]' }).textContent);

    init({
      pdbUrl: pdbUrl,
      ligandName: ligandName.trim(),
      interactions: interactions,
      hotspots: hotspots,
      lostResidues: lostResidues,
      gainedResidues: gainedResidues
    });
  }

  /* ---- Toolbar init ---- */
  function _initToolbar() {
    var resetBtn = document.getElementById('reset-view-btn');
    var pocketBtn = document.getElementById('pocket-view-btn');
    var hotspotBtn = document.getElementById('hotspot-view-btn');
    var surfaceBtn = document.getElementById('surface-toggle-btn');

    if (resetBtn) resetBtn.addEventListener('click', redraw);
    if (pocketBtn) pocketBtn.addEventListener('click', function () {
      if (viewer && _ligandName) { viewer.zoomTo({ resn: _ligandName }); viewer.render(); }
    });
    if (hotspotBtn) hotspotBtn.addEventListener('click', _focusHotspotsView);
    if (surfaceBtn) surfaceBtn.addEventListener('click', function () {
      surfaceOn = !surfaceOn;
      redraw();
    });
  }

  /* ---- Export ---- */
  window.PSCViewer = {
    init: init,
    initFromServer: initFromServer,
    redraw: redraw
  };

  // Auto-init toolbar on load
  document.addEventListener('DOMContentLoaded', _initToolbar);
})();
