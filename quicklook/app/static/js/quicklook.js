/* ============================================================================
 * QuickLook — shared client-side helpers
 *
 * Loaded by every GUI page (index, gallery, compare, tls_summary).
 * Owns:
 *  - Theme toggle + persistence (localStorage["ql-theme"])
 *  - Honours `prefers-color-scheme` on first visit
 *  - Wires the click handler on any element with id="themeToggle"
 *
 * Page-specific JS (form handling, WebSocket plumbing, table sorting,
 * etc.) stays in each template's inline <script>.
 * ============================================================================ */

(function () {
  'use strict';

  /* ---- Theme ----------------------------------------------------------- */
  function applyTheme(theme) {
    document.documentElement.setAttribute('data-theme', theme);
    const btn = document.getElementById('themeToggle');
    if (btn) btn.textContent = theme === 'dark' ? '☀' : '☾';
    try { localStorage.setItem('ql-theme', theme); } catch (e) {}
  }
  // Expose globally so page-specific scripts can call it if they need to.
  window.applyTheme = applyTheme;

  function initTheme() {
    let saved = null;
    try { saved = localStorage.getItem('ql-theme'); } catch (e) {}
    if (saved === 'light' || saved === 'dark') {
      applyTheme(saved);
      return;
    }
    if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
      applyTheme('dark');
    } else {
      // Force a light pass so the toggle glyph is correct from cold load.
      applyTheme('light');
    }
  }

  function wireThemeToggle() {
    const btn = document.getElementById('themeToggle');
    if (!btn) return;
    btn.addEventListener('click', () => {
      const cur = document.documentElement.getAttribute('data-theme');
      applyTheme(cur === 'dark' ? 'light' : 'dark');
    });
  }

  // Run as early as possible to minimise flash-of-unstyled-content.
  initTheme();
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', wireThemeToggle);
  } else {
    wireThemeToggle();
  }
})();
