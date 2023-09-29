import * as soda from "https://esm.run/@sodaviz/soda@0.13.1";
import * as d3 from "https://esm.run/d3-scale@4";

function run(data) {
  const LABEL_WIDTH = 500;
  let brushDomain = undefined;

  let timeoutTime = 50;
  let timeoutId = 0;

  let optionTimeoutTime = 100;
  let optionTimeoutId = 0;

  let traceColors = [
    "#4e79a7",
    "#f28e2c",
    "#e15759",
    "#76b7b2",
    "#59a14f",
    "#edc949",
    "#af7aa1",
    "#ff9da7",
    "#9c755f",
    "#bab0ab",
  ];

  // these maps contain the dom
  // nodes that control options
  let inputs = new Map();
  let values = new Map();

  let options = {
    colorScale: d3.scaleLinear().domain([0.0, 1.0]).range(["white", "black"]),
    segments: false,
    labels: true,
    confMin: 0.0,
    confMax: 1.0,
    confThresh: 0.0,
    aliThresh: 300,
    traceAtTop: true,
    filters: [],
    sidebarExpanded: false,
    sidebarMinWidth: 30,
    sidebarMaxWidth: 200,
  };

  let rowToQuery = [];
  let rowToLength = [];
  let rowToStrand = [];
  let rowToConf = [];

  let params = prepareData();
  let charts = initializeCharts();

  bindWidgets();

  render();

  function toggleSidebar() {
    let sidebar = document.querySelector("div#sidebar");
    if (options.sidebarExpanded === true) {
      sidebar.style.width = `${options.sidebarMinWidth}px`;
      options.sidebarExpanded = false;
    } else {
      sidebar.style.width = `${options.sidebarMaxWidth}px`;
      options.sidebarExpanded = true;
    }
  }

  function handleEvent(e) {
    let sliders = ["confMin", "confMax", "confThresh", "aliThresh"];
    let toggles = ["traceAtTop", "labels", "segments"];
    let text = ["filters"];
    let input = e.target.id;
    let value = e.target.value;

    if (sliders.indexOf(input) >= 0) {
      let valueNum = parseFloat(value);
      options[input] = valueNum;

      let valueNode = values.get(input);
      valueNode.textContent = value;
    } else if (toggles.indexOf(input) >= 0) {
      let toggleNode = inputs.get(input);
      let checked = toggleNode.checked;

      options[input] = checked;
    } else if (text.indexOf(input) >= 0) {
      options[input] = value.split("\n");
    } else {
      console.error("unknown input: ", input);
    }

    options.colorScale = d3
      .scaleLinear()
      .domain([options.confMin, options.confMax])
      .range(["white", "black"]);

    clearTimeout(optionTimeoutId);
    optionTimeoutId = window.setTimeout(() => {
      renderBottom(false);
    }, optionTimeoutTime);
  }

  function bindWidgets() {
    let toggle = document.querySelector("span#toggle");
    toggle.addEventListener("click", toggleSidebar);

    inputs.set("confMin", document.querySelector("input#confMin"));
    values.set("confMin", document.querySelector("div#confMin"));

    inputs.set("confMax", document.querySelector("input#confMax"));
    values.set("confMax", document.querySelector("div#confMax"));

    inputs.set("confThresh", document.querySelector("input#confThresh"));
    values.set("confThresh", document.querySelector("div#confThresh"));

    inputs.set("aliThresh", document.querySelector("input#aliThresh"));
    values.set("aliThresh", document.querySelector("div#aliThresh"));

    inputs.set("segments", document.querySelector("input#segments"));
    inputs.set("labels", document.querySelector("input#labels"));
    inputs.set("traceAtTop", document.querySelector("input#traceAtTop"));
    inputs.set("filters", document.querySelector("textArea#filters"));

    inputs.forEach((input) => input.addEventListener("input", handleEvent));
  }

  function initializeCharts() {
    // default chart settings that
    // are shared across all charts
    let chartConf = {
      selector: "div.viz",
      resizable: true,
      divOutline: "1px solid black",
      rowHeight: 16,
      padSize: 0,
    };

    let annChartConf = {
      ...chartConf,
      updateLayout(params) {
        let alignedWidthSort = (verts, graph) => {
          verts.sort((v1, v2) => {
            if (
              graph.getAnnotationFromId(v1).alignedWidth >
              graph.getAnnotationFromId(v2).alignedWidth
            ) {
              return -1;
            } else {
              return 1;
            }
          });
        };

        this.layout = soda.greedyGraphLayout(params.proxy, 0, alignedWidthSort);

        // sneaky: rewrite the layout object's row retrieval
        //         function so that it works for the annotations
        //         that the proxy annotations correspond to
        this.layout.row = function (d) {
          let id_tokens = d.a.id.split("-");
          let id = `${id_tokens[0]}-${id_tokens[1]}`;
          let row = this.rowMap.get(id);
          return row || 0;
        };
      },
      draw(params) {
        if (params.axis == true) {
          this.addAxis();
        }

        // fragments
        soda.chevronRectangle({
          chart: this,
          selector: "groups",
          annotations: params.aligned,
          orientation: (d) => d.a.strand,
          strokeColor: "black",
        });

        // joins;
        soda.line({
          chart: this,
          selector: "inner-join",
          annotations: params.inner,
        });

        // labels
        soda.dynamicText({
          chart: this,
          selector: "label",
          annotations: params.labels,
          fillColor: "black",
          y: (d) => this.layout.row(d) * this.rowHeight + 4,
          text: (d) => [d.a.label, "..."],
          textAnchor: "end",
        });

        soda.tooltip({
          chart: this,
          annotations: params.aligned,
          text: (d) => `${d.a.label}`,
        });
      },
    };

    let reference = new soda.Chart({
      upperPadSize: 25,
      ...annChartConf,
    });

    let aurora = new soda.Chart(annChartConf);
    let auroraZoom = new soda.Chart({
      ...annChartConf,
    });

    let genome = new soda.Chart({
      ...chartConf,
      draw(params) {
        soda.sequence({
          chart: this,
          selector: "genome",
          annotations: params.annotations,
        });
      },
    });

    let alignments = new soda.Chart({
      ...chartConf,
      zoomable: true,
      rowColors: ["whitesmoke", "white"],
      rowHeight: 30,

      updateLayout(params) {
        let queryIds = [...new Set(params.proxy.map((a) => a.queryId))];

        let traceQueryIds = [];

        if (options.traceAtTop) {
          traceQueryIds = [...new Set(params.trace.map((a) => a.queryId))];
        }

        let remainingQueryIds = queryIds.filter(
          (id) => traceQueryIds.indexOf(id) === -1,
        );

        let rowCount = 1;
        let dpRowToChartRow = new Map([[0, 0]]);

        let layoutFn = (id) => {
          // get the annotations with this id
          let rowAnn = params.proxy.filter((a) => a.queryId == id);

          // sub-group them by overlap
          let groups = soda.aggregateIntransitive({
            annotations: rowAnn,
            criterion: (a, b) => a.start < b.end && a.end > b.start,
          });

          // for each sub-group of overlapping annotations in the row
          let maxSubLayoutRowCount = 0;
          groups.forEach((g) => {
            let ann = g.annotations;
            let subLayout = soda.intervalGraphLayout(ann);
            maxSubLayoutRowCount = Math.max(
              maxSubLayoutRowCount,
              subLayout.rowCount,
            );
            // for each annotation in the sub-group
            ann.forEach((a) => {
              let subRow = subLayout.row({ a });
              dpRowToChartRow.set(a.row, rowCount + subRow);
            });
          });

          rowCount += maxSubLayoutRowCount;
        };

        traceQueryIds.forEach(layoutFn);
        remainingQueryIds.forEach(layoutFn);

        this.layout = {
          row: (d) => dpRowToChartRow.get(d.a.row),
          rowCount,
        };
      },
      draw(params) {
        this.clear();
        let domainWidth = this.domain[1] - this.domain[0];

        let domainFilter = (ann) =>
          ann.filter((a) => a.start < this.domain[1] && a.end > this.domain[0]);

        let y = (d) => this.rowHeight * this.layout.row(d) + 14;
        let x = (d) => this.xScale(d.a.start - 0.25);
        let width = (d) => this.xScale(d.a.end) - this.xScale(d.a.start - 0.5);

        if (options.segments) {
          // segments
          soda.rectangle({
            chart: this,
            selector: "segments",
            annotations: params.segments,
            row: 0,
            x,
            width,
            height: this.viewportHeightPx,
            fillColor: (d) => traceColors[d.a.traceIter],
            strokeColor: (d) => traceColors[d.a.traceIter],
            fillOpacity: 0.05,
            strokeWidth: 1,
          });
        } else {
          // note: this is a hack to place the segments
          //       <g> tag at the top until I can fix the
          //       bug in soda that doesn't clear <g> tags
          soda.rectangle({
            chart: this,
            selector: "segments",
            annotations: [],
          });
        }

        // proxy
        soda.rectangle({
          chart: this,
          selector: "proxy",
          annotations: params.proxy,
          y,
          width: (d) => this.xScale(d.a.end) - this.xScale(d.a.start + 1),
          height: 2,
        });

        // fragments
        soda.rectangle({
          chart: this,
          selector: "fragments",
          annotations: params.fragments,
          x,
          y,
          width,
          height: 12,
          fillColor: (d) => options.colorScale(d.a.conf),
          strokeColor: "black",
        });

        if (options.labels) {
          let exists = new Map();
          params.proxy.forEach((a) => exists.set(a.query, [false, a.row]));

          let annotations = domainFilter(params.proxy);

          annotations.forEach((a) => {
            if (this.xScale(a.start) < 200) {
              exists.set(a.query, [true, a.row]);
            }
          });

          let supplement = Array.from(exists.entries())
            .filter(([q, e]) => !e[0])
            .map(([q, e]) => {
              return { id: q, start: 1, end: 500, row: e[1], query: q };
            });

          annotations = annotations.concat(supplement);
          soda.dynamicText({
            chart: this,
            selector: "labels",
            annotations,
            x: (d) => Math.max(d.c.xScale(d.a.start), 0),
            text: (d) => [d.a.query, "..."],
          });
        }

        if (domainWidth < options.aliThresh) {
          let annotations = domainFilter(params.sequences);
          soda.sequence({
            chart: this,
            selector: "ali-seq",
            annotations,
            y,
            fillColor: (d) =>
              d.a.id[0] == "m" ? "green" : d.a.id[0] == "s" ? "orange" : "red",
          });
        }

        // trace
        soda.rectangle({
          chart: this,
          selector: "trace",
          annotations: params.trace,
          x,
          y,
          width,
          height: 12,
          fillColor: "none",
          strokeColor: (d) => traceColors[d.a.traceIter],
          strokeWidth: 3,
        });

        soda.tooltip({
          annotations: params.fragments,
          text: (d) =>
            `${d.a.query}: ` +
            `${d.a.queryStart.toLocaleString()}..${d.a.queryEnd.toLocaleString()}` +
            `<br>strand: ${d.a.strand}` +
            `<br>confidence: ${d.a.conf}`,
        });
      },

      postZoom() {
        clearTimeout(timeoutId);
        timeoutId = window.setTimeout(() => {
          this.draw({
            ...this.renderParams,
            updateDomain: false,
          });
        }, timeoutTime);
      },
    });

    alignments.render = function (params) {
      //this.resetTransform();

      let queryFilter = (a) => {
        let filters = options.filters;

        // don't ever filter the skip state
        if (a.row == 0) {
          return true;
        }

        let query = rowToQuery[a.row].toLowerCase();
        let conf = rowToConf[a.row];
        if (conf < options.confThresh) {
          return false;
        }

        let pass = true;
        filters.forEach((f) => {
          f = f.toLowerCase();
          if (f[0] !== "!") {
            if (!query.includes(f)) {
              pass = false;
              return;
            }
          } else {
            f = f.slice(1);
            if (query.includes(f)) {
              pass = false;
              return;
            }
          }
        });

        return pass;
      };

      let filteredParams = {
        start: params.start,
        end: params.end,
        segments: params.segments,
        fragments: params.fragments.filter(queryFilter),
        proxy: params.proxy.filter(queryFilter),
        sequences: params.sequences.filter(queryFilter),
        trace: params.trace.filter(queryFilter),
      };

      this.renderParams = filteredParams;
      this.updateLayout(filteredParams);
      this.updateRowCount(filteredParams);
      this.updateDimensions(filteredParams);
      if (params.updateDomain) {
        this.updateDomain(filteredParams);
      }
      this.draw(filteredParams);
      this.postRender(filteredParams);
    };
    return { reference, aurora, auroraZoom, genome, alignments };
  }

  function prepareAnn(ann) {
    let idCnt = 0;
    let proxy = [];
    let aligned = [];
    let inner = [];
    let left = [];
    let right = [];
    let labels = [];

    // for every group of joined annotations
    for (const group of ann) {
      // proxy -- for layout
      proxy.push({
        id: group.id,
        start: Math.min(group.visualStart, group.alignStart - LABEL_WIDTH),
        end: group.visualEnd,
        alignedWidth: group.alignEnd - group.alignStart + 1,
        label: group.query,
      });

      // aligned
      for (const record of group.aligned) {
        let tokens = record.split(",");
        aligned.push({
          id: tokens[0],
          start: parseInt(tokens[1]),
          end: parseInt(tokens[2]),
          strand: group.strand,
          label: group.query,
        });
      }

      // inner
      for (const record of group.inner) {
        let tokens = record.split(",");
        inner.push({
          id: tokens[0],
          start: parseInt(tokens[1]),
          end: parseInt(tokens[2]),
          queryLength: parseInt(tokens[3]),
        });
      }

      // left
      let tokens = group.left.split(",");
      left.push({
        id: tokens[0],
        start: parseInt(tokens[1]),
        end: parseInt(tokens[2]),
      });

      // right
      tokens = group.right.split(",");
      right.push({
        id: tokens[0],
        start: parseInt(tokens[1]),
        end: parseInt(tokens[2]),
      });

      // label
      labels.push({
        id: group.id,
        start: left[left.length - 1].end - 1,
        end: left[left.length - 1].end + LABEL_WIDTH,
        label: group.query,
      });
    }
    return {
      proxy,
      aligned,
      inner,
      left,
      right,
      labels,
    };
  }

  function prepareAli(ali) {
    let sequences = [];
    let proxy = [];
    let labelMap = new Map();
    let id = 0;
    let blank = "\u2000";
    for (const a of ali) {
      let tokens = a.split(",");
      let green = tokens[0].replace(/ /g, blank);
      let orange = tokens[1].replace(/ /g, blank);
      let start = parseInt(tokens[2]);
      let end = parseInt(tokens[3]);
      let query = tokens[4];
      let row = parseInt(tokens[5]);
      let queryId = parseInt(tokens[6]);
      let strand = tokens[7];

      labelMap.set(query, row);
      let common = { query, start, end, row, queryId, strand };

      sequences.push({
        id: `m-${id}`,
        ...common,
        sequence: green,
      });

      sequences.push({
        id: `s-${id}`,
        ...common,
        sequence: orange,
      });

      proxy.push({
        id: `ali-${id}`,
        ...common,
      });

      id++;
    }
    return { sequences, proxy };
  }

  function prepareTrace(traceStrings, targetStart) {
    let trace = [];
    for (const [iter, line] of traceStrings.entries()) {
      let iterStrings = line.split("|");
      for (const [idx, seg] of iterStrings.entries()) {
        let tokens = seg.split(",");
        let start = parseInt(tokens[0]) + targetStart;
        let end = parseInt(tokens[1]) + targetStart;
        let queryId = parseInt(tokens[2]);
        let row = parseInt(tokens[3]);
        let query = tokens[4];
        trace.push({
          id: `trace-${iter}-${idx}`,
          traceIter: iter,
          start,
          end,
          queryId,
          row,
          query,
        });
      }
    }

    return { trace };
  }

  function prepareFragments(fragmentStrings, targetStart) {
    let fragments = [];
    for (const [iter, line] of fragmentStrings.entries()) {
      let iterStrings = line.split("|");
      for (const [idx, seg] of iterStrings.entries()) {
        let tokens = seg.split(",");
        let start = parseInt(tokens[0]) + targetStart;
        let end = parseInt(tokens[1]) + targetStart;
        let row = parseInt(tokens[2]);
        let conf = parseFloat(tokens[3]);
        let query = tokens[4];
        let queryStart = parseInt(tokens[5]);
        let queryEnd = parseInt(tokens[6]);
        let strand = tokens[7];
        fragments.push({
          id: `fragment-${iter}-${idx}`,
          traceIter: iter,
          start,
          end,
          conf,
          row,
          query,
          queryStart,
          queryEnd,
          strand,
        });
      }
    }

    return { fragments };
  }

  function prepareSegments(segmentStrings, targetStart) {
    let segments = [];
    for (const [iter, line] of segmentStrings.entries()) {
      let iterStrings = line.split("|");
      for (const [idx, seg] of iterStrings.entries()) {
        let tokens = seg.split(",");
        let start = parseInt(tokens[0]) + targetStart;
        let end = parseInt(tokens[1]) + targetStart;
        segments.push({
          id: `segment-${iter}-${idx}`,
          traceIter: iter,
          start,
          end,
        });
      }
    }

    return { segments };
  }

  function prepareData() {
    let coords = {
      start: data.targetStart - LABEL_WIDTH,
      end: data.targetEnd,
    };

    let aurora = {
      ...coords,
      ...prepareAnn(data.auroraAnn),
    };

    let reference = {
      axis: true,
      ...coords,
      ...prepareAnn(data.referenceAnn),
    };

    let genome = {
      ...coords,
      annotations: [
        {
          id: "genome-sequence",
          start: data.targetStart,
          end: data.targetEnd + 1,
          sequence: data.targetSeq,
        },
      ],
    };

    let alignments = {
      ...coords,
      ...prepareAli(data.alignmentStrings),
      ...prepareTrace(data.traceStrings, data.targetStart),
      ...prepareFragments(data.fragmentStrings, data.targetStart),
      ...prepareSegments(data.segmentStrings, data.targetStart),
    };

    rowToLength[0] = data.targetEnd - data.targetStart + 1;
    alignments.proxy.forEach((a) => {
      rowToQuery[a.row] = a.query;
      rowToLength[a.row] = a.end - a.start + 1;
      rowToStrand[a.row] = a.strand;
    });

    for (let row = 0; row < rowToQuery.length; row++) {
      let fragments = alignments.fragments.filter((f) => f.row == row);
      let confTotal = fragments.reduce(
        (acc, f) => acc + f.conf * (f.end - f.start + 1),
        0,
      );
      rowToConf[row] = confTotal / rowToLength[row];
    }

    return { aurora, reference, genome, alignments };
  }

  function initializeBrush() {
    let chart = charts["aurora"];
    chart.viewportSelection.call(
      soda.internalD3
        .brushX()
        .extent([
          [0, 0],
          [chart.viewportWidthPx, chart.viewportHeightPx + 1],
        ])
        .on("start", () => {})
        .on("brush", () => {
          let brushRange = soda.internalD3.event.selection;
          brushDomain = [
            Math.round(chart.xScale.invert(brushRange[0])),
            Math.round(chart.xScale.invert(brushRange[1])),
          ];
        })
        .on("end", () => renderBottom()),
    );
  }

  function renderBottom(updateDomain = true) {
    if (brushDomain == undefined) {
      return;
    }

    let coords = {
      start: brushDomain[0],
      end: brushDomain[1],
    };

    charts.genome.render({
      ...params.genome,
      ...coords,
    });

    charts.auroraZoom.render({
      ...params.aurora,
      ...coords,
    });

    charts.alignments.render({
      updateDomain,
      ...params.alignments,
      ...coords,
    });
  }

  function render() {
    charts.reference.render(params.reference);
    charts.aurora.render(params.aurora);
    initializeBrush();

    let zoomSync = new soda.ZoomSyncer();
    zoomSync.add([charts.auroraZoom, charts.genome, charts.alignments]);
  }
}
