/* ============================================================
   api.js — All fetch/API calls, unified error handling
   ============================================================ */

(function () {
  'use strict';

  var FETCH_TIMEOUT = 60000; // 60s

  function _fetch(url, options) {
    var controller = new AbortController();
    var timer = setTimeout(function () { controller.abort(); }, FETCH_TIMEOUT);
    var opts = Object.assign({ signal: controller.signal }, options || {});
    return fetch(url, opts).then(function (r) {
      clearTimeout(timer);
      if (!r.ok) throw new Error('HTTP ' + r.status);
      return r.json();
    }).catch(function (err) {
      clearTimeout(timer);
      if (err.name === 'AbortError') throw new Error('Request timed out');
      throw err;
    });
  }

  /* ---- RCSB Fetch ---- */
  function fetchPDB(pdbId) {
    var formData = new FormData();
    formData.append('pdb_id', pdbId.toUpperCase());
    return _fetch('/fetch_pdb', { method: 'POST', body: formData });
  }

  /* ---- Analyze (AJAX) ---- */
  function analyze(formData) {
    formData.append('skipLigand', AppState.get('ui').skipLigand ? 'true' : 'false');
    return _fetch('/analyze_async', {
      method: 'POST',
      body: formData,
      headers: { 'X-Requested-With': 'XMLHttpRequest' }
    });
  }

  /* ---- Mutation Scan ---- */
  function mutationScan(formData) {
    return _fetch('/mutation_scan', {
      method: 'POST',
      body: formData,
      headers: { 'X-Requested-With': 'XMLHttpRequest' }
    });
  }

  /* ---- Comparison ---- */
  function compare(formData) {
    return _fetch('/compare', {
      method: 'POST',
      body: formData,
      headers: { 'X-Requested-With': 'XMLHttpRequest' }
    });
  }

  /* ---- Stats ---- */
  function fetchStats() {
    return _fetch('/api/stats');
  }

  function fetchRecent() {
    return _fetch('/api/recent_analyses');
  }

  /* ---- Export ---- */
  window.API = {
    fetchPDB: fetchPDB,
    analyze: analyze,
    mutationScan: mutationScan,
    compare: compare,
    fetchStats: fetchStats,
    fetchRecent: fetchRecent
  };
})();
