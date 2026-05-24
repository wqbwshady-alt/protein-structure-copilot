/* ============================================================
   state.js — Single source of truth for all application state

   Usage:
     AppState.set('analysis', { ready: true, name: '1HSG.pdb' })
     AppState.on('analysis', function(data) { ... })
     AppState.get('analysis')  // { ready: true, name: '1HSG.pdb' }
   ============================================================ */

(function () {
  'use strict';

  var _state = {};
  var _listeners = {};

  /* ---- Empty structure template ---- */
  function emptyStructure() {
    return {
      name: null, size: null, source: null,
      fileObject: null, pdbFilename: null, pdbText: null, pdbId: null,
      ready: false
    };
  }

  /* ---- Initial state ---- */
  function _initialState() {
    return {
      mode: 'single',
      analysis: emptyStructure(),
      wt: emptyStructure(),
      mutant: emptyStructure(),
      ligandDetection: { status: 'unknown', ligands: [], primaryLigand: '' },
      results: {
        single: null,
        mutation: null,
        compare: null
      },
      viewer: {
        ligandName: '',
        interactions: [],
        hotspots: [],
        lostResidues: [],
        gainedResidues: [],
        pdbUrl: null
      },
      stats: { totalAnalyses: 0, recent: [] },
      ui: { loading: false, error: null }
    };
  }

  _state = _initialState();

  /* ---- Public API ---- */

  /** Get a top-level key or a dotted path */
  function get(path) {
    if (!path) return _state;
    var keys = path.split('.');
    var value = _state;
    for (var i = 0; i < keys.length; i++) {
      if (value == null) return undefined;
      value = value[keys[i]];
    }
    return value;
  }

  /** Set a top-level key with deep merge */
  function set(key, value) {
    var old = _state[key];
    if (typeof value === 'object' && value !== null && !Array.isArray(value)) {
      _state[key] = Object.assign({}, _state[key], value);
    } else {
      _state[key] = value;
    }
    _notify(key, _state[key], old);
    return _state[key];
  }

  /** Replace entire state at key */
  function replace(key, value) {
    var old = _state[key];
    _state[key] = value;
    _notify(key, value, old);
    return _state[key];
  }

  /** Subscribe to changes on a top-level key */
  function on(key, callback) {
    if (!_listeners[key]) _listeners[key] = [];
    _listeners[key].push(callback);
    // Return unsubscribe function
    return function () {
      _listeners[key] = (_listeners[key] || []).filter(function (cb) { return cb !== callback; });
    };
  }

  /** Reset to initial state */
  function reset() {
    _state = _initialState();
  }

  function _notify(key, newVal, oldVal) {
    (_listeners[key] || []).forEach(function (cb) {
      try { cb(newVal, oldVal); } catch (e) { console.error('State listener error:', e); }
    });
    // Also notify 'all' listeners
    (_listeners['*'] || []).forEach(function (cb) {
      try { cb(key, newVal, oldVal); } catch (e) { console.error('State listener error:', e); }
    });
  }

  /* ---- Export ---- */
  window.AppState = {
    get: get,
    set: set,
    replace: replace,
    on: on,
    reset: reset,
    emptyStructure: emptyStructure
  };
})();
