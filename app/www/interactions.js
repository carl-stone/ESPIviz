(function () {
  window.ESPIviz = window.ESPIviz || {};
  window.ESPIviz.resizeDEPlot = function (plot) {
    if (plot._espivizResizeObserver) plot._espivizResizeObserver.disconnect();
    var ticks = plot.layout.xaxis.tickvals;
    var lastWidth = 0;
    function resize() {
      var width = plot.clientWidth;
      if (!width || width === lastWidth) return;
      lastWidth = width;
      var compact = width < 850;
      var narrow = width < 480;
      var left = narrow ? 76 : 86;
      var right = 16;
      var layout = {
        'font.size': narrow ? 12 : 13,
        'xaxis.title.font.size': narrow ? 12 : 15,
        'yaxis.title.font.size': narrow ? 12 : 15,
        'xaxis.tickfont.size': narrow ? 11 : 13,
        'yaxis.tickfont.size': narrow ? 11 : 13,
        'xaxis.tickangle': 0,
        'xaxis.nticks': narrow ? 4 : 7,
        'legend.orientation': compact ? 'v' : 'h',
        'legend.font.size': narrow ? 12 : 13,
        'legend.title.font.size': narrow ? 12 : 13,
        'legend.x': compact ? (8 - left) / Math.max(1, width - left - right) : 0,
        'legend.y': compact ? -0.18 : -0.2,
        'margin': {l: left, r: right, t: 18, b: compact ? 178 : 124}
      };
      if (ticks) {
        var visibleTicks = ticks.filter(function (_, i) { return !narrow || i % 2 === 0; });
        layout['xaxis.tickvals'] = visibleTicks;
        layout['xaxis.ticktext'] = visibleTicks.map(function (value) {
          if (!narrow) return value.toLocaleString('en-US');
          if (value >= 1000000) return value / 1000000 + 'M';
          if (value >= 1000) return value / 1000 + 'k';
          return String(value);
        });
      }
      window.Plotly.relayout(plot, layout);
    }
    plot._espivizResizeObserver = new ResizeObserver(resize);
    plot._espivizResizeObserver.observe(plot);
    resize();
  };
})();

(function () {
  function closeMenus(except) {
    document.querySelectorAll('.tool-menu[open]').forEach(function (menu) {
      if (menu !== except) menu.open = false;
    });
  }
  document.addEventListener('toggle', function (event) {
    var details = event.target;
    if (!(details instanceof HTMLDetailsElement)) return;
    if (details.matches('.tool-menu') && details.open) closeMenus(details);
    if (details.matches('.data-disclosure')) {
      window.dispatchEvent(new Event('resize'));
    }
  }, true);
  document.addEventListener('click', function (event) {
    if (!event.target.closest('.tool-menu')) closeMenus(null);
  });
  document.addEventListener('keydown', function (event) {
    if (event.key !== 'Escape') return;
    var openMenu = document.querySelector('.tool-menu[open]');
    if (openMenu) {
      closeMenus(null);
      openMenu.querySelector('summary').focus();
    }
  });
  document.addEventListener('bslib.card', function (event) {
    var card = event.target;
    var button = card.querySelector('.bslib-full-screen-enter');
    if (!button) return;
    button.setAttribute('aria-expanded', String(event.detail.fullScreen));
    if (!event.detail.fullScreen) button.focus();
  });
})();

(function () {
  function whenShiny() {
    if (!window.Shiny) return setTimeout(whenShiny, 100);
    Shiny.addCustomMessageHandler('control-disabled', function (data) {
      var control = document.getElementById(data.id);
      if (control) { control.disabled = data.disabled; control.setAttribute('aria-disabled', String(data.disabled)); if (control.tagName === 'A') { control.classList.toggle('disabled', data.disabled); control.tabIndex = data.disabled ? -1 : 0; } }
    });
    Shiny.addCustomMessageHandler('download-label', function (data) {
      var control = document.getElementById(data.id);
      if (control) { var icon = control.querySelector('i'); control.textContent = ' ' + data.label; if (icon) control.prepend(icon); }
    });
    Shiny.addCustomMessageHandler('copy-view', async function (data) {
      var status = document.getElementById(data.statusId);
      var url = window.location.origin + window.location.pathname + data.query;
      try {
        if (navigator.clipboard && window.isSecureContext) await navigator.clipboard.writeText(url);
        else {
          var field = document.createElement('textarea'); field.value = url;
          field.style.position = 'fixed'; field.style.opacity = '0'; document.body.appendChild(field);
          field.select(); var ok = document.execCommand('copy'); field.remove();
          if (!ok) throw new Error('copy');
        }
        if (status) status.textContent = 'View link copied. It includes this data version and all supported view settings.';
      } catch (_) {
        if (status) { status.textContent = 'Copy this view link: '; var field = document.createElement('input'); field.value = url; field.readOnly = true; field.setAttribute('aria-label', 'View link'); status.appendChild(field); field.focus(); field.select(); }
      }
    });
  }
  whenShiny();
  function improveTables() {
    document.querySelectorAll('i[role="presentation"][aria-label]').forEach(function (icon) { icon.removeAttribute("aria-label"); icon.setAttribute("aria-hidden", "true"); });
    document.querySelectorAll('.dataTables_wrapper').forEach(function (wrapper) {
      // With scrollX, DataTables puts a header-only clone before the body table.
      var table = wrapper.querySelector('.dataTables_scrollBody table.dataTable') || wrapper.querySelector('table.dataTable');
      if (!table) return;
      var heading = wrapper.closest('.card')?.querySelector('h3, .card-header');
      var name = heading ? heading.textContent.trim() : 'Result';
      table.setAttribute('aria-label', name + ' table');
      wrapper.querySelectorAll('.dataTables_filter input').forEach(function (input) { input.setAttribute('aria-label', 'Search ' + name.toLowerCase() + ' table'); });
      wrapper.querySelectorAll('thead tr').forEach(function (row) {
        row.querySelectorAll('input,select').forEach(function (input, index) {
          var headers = table.querySelectorAll('thead tr:first-child th');
          input.setAttribute('aria-label', 'Filter ' + (headers[index]?.textContent.trim() || name));
        });
      });
      if (table.dataset.rowSelectable === 'true') {
        table.querySelectorAll('tbody tr').forEach(function (row) {
          if (row.querySelector('.dataTables_empty')) return;
          row.tabIndex = 0;
          row.setAttribute('aria-label', 'Result ' + row.cells[0]?.textContent.trim() + '. Press Enter to select.');
          if (!row.dataset.keyboardBound) {
            row.dataset.keyboardBound = 'true';
            row.addEventListener('keydown', function (event) {
              if (event.key !== 'Enter' && event.key !== ' ') return;
              event.preventDefault();
              // DT's row selection is bound to mousedown, including server-side indices.
              row.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
            });
          }
        });
      }
    });
  }
  if (window.jQuery) { window.jQuery(document).on('draw.dt shiny:value', function () { setTimeout(improveTables, 20); }); window.jQuery(document).on('shown.bs.tab', '#main_nav a', function () { window.scrollTo({top: 0, behavior: 'instant'}); }); }
  document.addEventListener('DOMContentLoaded', improveTables);
  document.addEventListener('toggle', function (event) {
    if (event.target.matches && event.target.matches('.data-disclosure')) setTimeout(improveTables, 20);
  }, true);
})();
