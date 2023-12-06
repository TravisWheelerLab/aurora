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
    sidebarMinWidth: 30,
    sidebarMaxWidth: 200,
  };

  let state = {
    sidebarExpanded: false,
    colorScale: d3.scaleLinear().domain([0.0, 1.0]).range(["white", "black"]),
    labels: true,
    traceAtTop: true,
    assemblyAtTop: false,
    confThresh: 0.0,
    confMin: 0.0,
    confMax: 1.0,
    aliThresh: 300,
    regex: undefined,
    traceIteration: undefined,
    numTraceIterations: undefined,
    onlyTrace: true,
    traceRowsByIter: undefined,
  }

  let rowToQuery = [];

  let params = prepareData();
  let charts = initializeCharts();

  bindWidgets();

  render();

  function toggleSidebar() {
    let sidebar = document.querySelector("div#sidebar");
    if (state.sidebarExpanded === true) {
      sidebar.style.width = `${options.sidebarMinWidth}px`;
      state.sidebarExpanded = false;
    } else {
      sidebar.style.width = `${options.sidebarMaxWidth}px`;
      state.sidebarExpanded = true;
    }
  }

  function handleEvent(e) {
    let toggles = ["traceAtTop", "labels", "onlyTrace"];
    let numeric = ["confThresh", "aliThresh"]
    let text = ["regex"];
    let traceButtons = [];
    for (let i = 0; i < state.numTraceIterations; i++) {
      traceButtons.push(`traceButton${i}`);
    }

    let input = e.target.id;
    let value = e.target.value;

    if (toggles.indexOf(input) >= 0) {
      let toggleNode = inputs.get(input);
      let checked = toggleNode.checked;
      state[input] = checked;
    } else if (numeric.indexOf(input) >= 0) {
      state[input] = parseFloat(value);
    } else if (text.indexOf(input) >= 0) {
      if (value != "") {
        state[input] = new RegExp(value);
      } else {
        state[input] = undefined;
      }
    } else if (traceButtons.indexOf(input) >= 0) {
      let traceIter = traceButtons.indexOf(input);
      state.traceIteration = traceIter;
    } else {
      console.error("unknown input: ", input);
    }

    state.colorScale = d3
      .scaleLinear()
      .domain([state.confMin, state.confMax])
      .range(["white", "black"]);

    clearTimeout(optionTimeoutId);
    optionTimeoutId = window.setTimeout(() => {
      renderBottom(false);
    }, optionTimeoutTime);
  }

  function bindWidgets() {
    let toggle = document.querySelector("span#toggle");
    toggle.addEventListener("click", toggleSidebar);

    let confThresh = document.querySelector("input#confThresh");
    confThresh.value = state.confThresh;
    inputs.set("confThresh", confThresh);

    let aliThresh = document.querySelector("input#aliThresh");
    aliThresh.value = state.aliThresh;
    inputs.set("aliThresh", aliThresh);

    inputs.set("onlyTrace", document.querySelector("input#onlyTrace"));
    inputs.set("labels", document.querySelector("input#labels"));
    inputs.set("traceAtTop", document.querySelector("input#traceAtTop"));
    inputs.set("regex", document.querySelector("input#regex"));

    // grab the entire sidebar
    let innerSidebarDiv = document.querySelector("div.inner-sidebar");

    // create a div for the trace selection widget
    let traceSelectionWidget = document.createElement("div");
    traceSelectionWidget.className = "widget-vertical";

    // shove that widget into the sidebar
    innerSidebarDiv.appendChild(traceSelectionWidget);

    // programmatically create a radio button for each trace iteration
    for (let i = 0; i < state.numTraceIterations; i++) {
      let radioDiv = document.createElement('div');
      let radioButton = document.createElement('input');
      radioButton.type = 'radio';
      radioButton.name = 'traceButtonGroup';
      radioButton.id = `traceButton${i}`;

      let label = document.createElement('label');
      label.htmlFor = `traceButton${i}`;
      label.appendChild(document.createTextNode(`Trace ${i}`));

      radioDiv.appendChild(radioButton);
      radioDiv.appendChild(label);

      traceSelectionWidget.append(radioDiv);
    }
    document.querySelector(`#traceButton0`).checked = true;

    inputs.set("trace", traceSelectionWidget);

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
        this.layout.row = function(d) {
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
      upperPadSize: 25,
      draw(params) {
        this.addAxis();

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

        if (state.traceAtTop) {
          let allTrace = params.conclusiveTrace.concat(params.ambiguousTrace);
          allTrace.sort((a, b) => a.start - b.start);
          traceQueryIds = [...new Set(allTrace.map((a) => a.queryId))];
        }

        let remainingQueryIds = queryIds.filter(
          (id) => traceQueryIds.indexOf(id) === -1,
        ).filter(
          (id) => assemblyQueryIds.indexOf(id) === -1,
        );

        // the first two rows are reserved for the 
        // skip state and the tandem repeat state
        let rowCount = 2;
        let dpRowToChartRow = new Map([[0, 0]]);
        params.tandemRepeats.forEach((r) => {
          dpRowToChartRow.set(r.row, 1);
        }) 
      
        let layoutFn = (id) => {
          let queryAssemblies = params.assemblies.filter((a) => a.queryId == id);
          let subLayout = soda.intervalGraphLayout(queryAssemblies);

          queryAssemblies.forEach((a) => {
            let subRow = subLayout.row({ a });
            dpRowToChartRow.set(a.row, rowCount + subRow)
          });

          rowCount += subLayout.rowCount;
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

        // assemblies
        soda.rectangle({
          chart: this,
          selector: "assembly",
          annotations: params.assemblies.filter((a) => a.size > 1),
          y,
          width: (d) => this.xScale(d.a.end) - this.xScale(d.a.start + 1),
          height: 2,
          fillColor: (d) => {
            if (params.competedAssemblyRows[state.traceIteration].indexOf(d.a.row) >= 0) {
              return "red";
            }
            else if (params.unresolvedAssemblyRows[state.traceIteration].indexOf(d.a.row) >= 0) {
              return "orange";
            }
            else if (params.resolvedAssemblyRows[state.traceIteration].indexOf(d.a.row) >= 0) {
              return "green";
            }
            else {
              return "black";
            }
          },
        });

        // alignments
        soda.rectangle({
          chart: this,
          selector: "alignments",
          annotations: params.sequences,
          x,
          y,
          width,
          height: 12,
          fillOpacity: 0.05,
          strokeColor: "black",
        });

        // tandem repeats
        soda.rectangle({
          chart: this,
          selector: "tandem-repeats",
          annotations: params.tandemRepeats,
          x,
          y,
          width,
          height: 12,
          fillOpacity: 0.05,
          strokeColor: "black",
        });

        if (state.labels) {
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

        if (domainWidth < state.aliThresh) {
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
          selector: "conclusive-trace",
          annotations: params.conclusiveTrace,
          x,
          y,
          width,
          height: 3,
          fillColor: "none",
          strokeColor: "blue",
          strokeWidth: 3,
        });

        soda.rectangle({
          chart: this,
          selector: "ambiguous-trace",
          annotations: params.ambiguousTrace,
          x,
          y,
          width,
          height: 3,
          fillColor: "none",
          strokeColor: "orange",
          strokeWidth: 3,
        });

        // soda.tooltip({
        //   annotations: params.sequences,
        //   text: (d) =>
        //     `${d.a.query}: ` +
        //     `${d.a.queryStart.toLocaleString()}..${d.a.queryEnd.toLocaleString()}` +
        //     `<br>chrom: ${d.a.start.toLocaleString()}..${d.a.end.toLocaleString()}` +
        //     `<br>strand: ${d.a.strand}` +
        //     `<br>confidence: ${d.a.conf}`,
        // });
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

    alignments.render = function(params) {
      //this.resetTransform();

      let queryFilter = (a) => {
        // don't ever filter the skip state
        if (a.row == 0) {
          return true;
        }

        if (a.row > rowToQuery.length - 1) {
          return true;
        }

        if (state.onlyTrace) {
          if (state.traceRowsByIter[state.traceIteration].indexOf(a.row) < 0) {
            return false;
          }
        }
        

        let query = rowToQuery[a.row].toLowerCase();

        if (state.regex != undefined)  {
          return state.regex.test(query);
        }

        return true; 
      };

      let filteredParams = {
        ...params,
        assemblies: params.assemblies.filter(queryFilter),
        proxy: params.proxy.filter(queryFilter),
        sequences: params.sequences.filter(queryFilter),
        ambiguousTrace: params.ambiguousTrace[state.traceIteration].filter(queryFilter),
        conclusiveTrace: params.conclusiveTrace[state.traceIteration].filter(queryFilter),
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
    let traces = [];
    for (const [iter, line] of traceStrings.entries()) {
      let iterStrings = line.split("|");
      let iterTrace = [];
      for (const [idx, seg] of iterStrings.entries()) {
        if (seg == "") {
          continue;
        }
        let tokens = seg.split(",");
        let start = parseInt(tokens[0]) + targetStart;
        let end = parseInt(tokens[1]) + targetStart;
        let queryId = parseInt(tokens[2]);
        let row = parseInt(tokens[3]);
        iterTrace.push({
          id: `trace-${iter}-${idx}`,
          traceIter: iter,
          start,
          end,
          queryId,
          row,
        });
      }
      traces.push(iterTrace);
    }
    return traces;
  }



  function prepareAssemblies(assemblyStrings) {
    let assemblies = [];
    for (const [idx, seg] of assemblyStrings.entries()) {
      let tokens = seg.split(",");
      let start = parseInt(tokens[0]);
      let end = parseInt(tokens[1]);
      let queryId = parseInt(tokens[2]);
      let size = parseInt(tokens[3]);
      let row = parseInt(tokens[4]);
      assemblies.push({
        id: `assembly-${idx + 1}`,
        queryId,
        start,
        end,
        size,
        row,
      });
    }
    return { assemblies };
  }

  function prepareTandemRepeats(tandemRepeatStrings) {
    let tandemRepeats = [];
    for (const seg of tandemRepeatStrings) {
      let tokens = seg.split(",");
      let start = parseInt(tokens[0]);
      let end = parseInt(tokens[1]);
      let consensus = tokens[2];
      let row = parseInt(tokens[3]);
      tandemRepeats.push({
        id: `tr-${row}`,
        start,
        end,
        consensus,
        row,
      });
    }
    return { tandemRepeats };
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
      ...prepareAssemblies(data.assemblyStrings),
      ...prepareTandemRepeats(data.tandemRepeatStrings),
      ambiguousTrace: prepareTrace(data.ambiguousTraceStrings, data.targetStart),
      conclusiveTrace: prepareTrace(data.conclusiveTraceStrings, data.targetStart),
      resolvedAssemblyRows: data.resolvedAssemblyRows,
      unresolvedAssemblyRows: data.unresolvedAssemblyRows,
      competedAssemblyRows: data.competedAssemblyRows,
    };

    alignments.proxy.forEach((a) => rowToQuery[a.row] = a.query);

    state.traceIteration = 0;
    state.numTraceIterations = alignments.ambiguousTrace.length;

    let traceRowsByIter = [];
    for (let i = 0; i < state.numTraceIterations; i++) {
      let rows = [];
      alignments.ambiguousTrace[i].forEach((a) => rows.push(a.row));
      alignments.conclusiveTrace[i].forEach((a) => rows.push(a.row));

      traceRowsByIter.push([...new Set(rows)]);
    }
    state.traceRowsByIter = traceRowsByIter;

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
        .on("start", () => { })
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
