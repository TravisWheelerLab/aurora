SODA_TARGET;
FZSTD_TARGET;

let HTML_TEMPLATE = `
  HTML_TARGET
`;

async function load_b64_gzip_json(text) {
  let rawText;
  if (Uint8Array.fromBase64) {
    rawText = Uint8Array.fromBase64(text);
  } else {
    text = btoa(text);
    rawText = new Uint8Array(text.length);
    for (let i = 0; i < text.length; i++)
      rawText[i] = text.charCodeAt(i);
  }

  const json = new TextDecoder().decode(fzstd.decompress(rawText));
  return JSON.parse(json);
}

async function bootstrap() {
  let confidences = document.getElementById("confidences");
  let data_elm = document.getElementById("data");

  if (!data_elm) {
    console.error("Unable to find annotation data.")
    document.body.innerText = "Unable to find annotation data.";
    return;
  }

  let data = await load_b64_gzip_json(data_elm.textContent);

  if(confidences) {
    data.alignmentConfidences = await load_b64_gzip_json(confidences.textContent);
  }

  window.auroraData = data;

  let parser = new DOMParser();
  let relative_path = document.getElementById("code").src.split("/").slice(0, -1);
  let html_template = HTML_TEMPLATE
    .replace("ICON_PATH_TARGET", [...relative_path, "icon.svg"].join("/"))
    .replace("INDEX_PATH_TARGET", [...relative_path, "index.html"].join("/"));
  let template_doc = parser.parseFromString(html_template, "text/html");

  // Update the page...
  document.head.append(...template_doc.head.children);
  document.body = template_doc.body;

  run(data);
}


function run(data) {
  document
    .querySelector(".container")
    .addEventListener("wheel", function (event) {
      if (event.ctrlKey) {
        event.preventDefault();
      }
    });

  const LABEL_WIDTH = 500;
  let brushDomain = undefined;

  let timeoutTime = 50;
  let timeoutId = 0;

  let optionTimeoutTime = 100;
  let optionTimeoutId = 0;

  let classColorCache = {
    "sine": "#1f77b4",
    "line": "#ff7f0e",
    "ltr": "#2ca02c",
    "dna": "#d62728",
    "simple": "#9467bd",
    "low_complexity": "#8c564b",
    "satellite": "#e377c2",
    "rna": "#7f7f7f",
    "other": "#bcbd22",
    "unknown": "#17becf",
    "unspecified": "#17becf",
    "simple_repeat": "#d13d8e",
    "tandem-repeat": "#d13d8e",
  }

  function stringToColor(s, heavyChars = 5) {
      let validChars = 69;
      let hash = 0;
      let revHash = 0;
      let multiplier = 360 * (2 ** heavyChars);
      let stepDivisor = 2;
      let range = 360;

      let m = multiplier / (validChars - 1);

      for(let i = 0; i < s.length; i++) {
          hash = (hash + charToNumber(s.charAt(i)) * m) % range;
          revHash = (revHash + charToNumber(s.charAt(s.length - (1 + i))) * m) % range;
          m /= stepDivisor;
      }

      let sat = 60 + (revHash % 30);
      let light = 30 + ((revHash * 30) % 50);

      return "hsl(" + hash + ", " + sat + "%, " + light + "%)";
  }


  function charToNumber(char) {
      if("a" <= char && char <= "z") {
          return (char.charCodeAt(0) - "a".charCodeAt(0)) + 1;
      }
      else {
          return 0;
      }
  }

  function classColorFromName(name) {
    for (let splitChar of ["#", "_"]) {
      let nameSplit = name.split(splitChar);

      if (nameSplit.length > 1) {
        let c = nameSplit[1].split("/")[0].toLowerCase();
        if(!(c in classColorCache)) {
          classColorCache[c] = stringToColor(c);
        }
        return classColorCache[c];
      }
    }
    return classColorCache["unknown"];
  }

  // this maps contain the dom
  // nodes that control options
  let inputs = new Map();

  let options = {
    sidebarMinWidth: 30,
    sidebarMaxWidth: 200,
    conclusiveColor: "#4e79a7",
    ambiguousColor: "#f28e2c",
    resolvedColor: "green",
    unresolvedColor: "orange",
    competedColor: "red",
    confidenceSegmentColor: "purple",
  };

  let state = {
    sidebarExpanded: false,
    labels: true,
    traceAtTop: true,
    onlySelected: true,
    assemblyAtTop: false,
    confThresh: 0.0,
    aliThresh: 300,
    regex: undefined,
    traceIteration: undefined,
    numTraceIterations: undefined,
    onlyTrace: true,
    showSegments: false,
    showInversions: false,
    showInactive: false,
    traceRowsByIter: undefined,
  };

  let rowToQuery = [];
  let rowToConf = [];

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
    let toggles = [
      "traceAtTop",
      "onlySelected",
      "labels",
      "onlyTrace",
      "showInactive",
      "showSegments",
      "showInversions",
    ];
    let numeric = ["confThresh", "aliThresh"];
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

    clearTimeout(optionTimeoutId);
    optionTimeoutId = window.setTimeout(() => {
      if (input == "showInversions") {
        charts.reference.render(params.reference);
        charts.aurora.render(params.aurora);
      }
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

    inputs.set("labels", document.querySelector("input#labels"));
    inputs.set("traceAtTop", document.querySelector("input#traceAtTop"));
    inputs.set("onlySelected", document.querySelector("input#onlySelected"));
    inputs.set("onlyTrace", document.querySelector("input#onlyTrace"));
    inputs.set("showInactive", document.querySelector("input#showInactive"));
    inputs.set("showSegments", document.querySelector("input#showSegments"));
    inputs.set(
      "showInversions",
      document.querySelector("input#showInversions"),
    );
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
      let radioDiv = document.createElement("div");
      let radioButton = document.createElement("input");
      radioButton.type = "radio";
      radioButton.name = "traceButtonGroup";
      radioButton.id = `traceButton${i}`;

      let label = document.createElement("label");
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

        function classColor(d) {
          return classColorFromName(d.a.label);
        }

        soda.rectangle({
          chart: this,
          selector: "inversion-highlight",
          annotations: state.showInversions ? params.inversions : [],
          fillColor: "yellow",
          fillOpacity: 0.5,
          y: 0,
          height: this.viewportHeightPx,
        });

        // fragments
        soda.chevronRectangle({
          chart: this,
          selector: "groups",
          annotations: params.aligned,
          orientation: (d) => d.a.strand,
          strokeColor: classColor,
          fillColor: (d) => {
            return state.showInversions && d.a.in_inversion
              ? "yellow"
              : "white";
          },
          strokeWidth: 2,
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

    let referenceZoom = new soda.Chart(annChartConf);
    let auroraZoom = new soda.Chart(annChartConf);

    let genome = new soda.Chart({
      ...chartConf,
      upperPadSize: 25,
      updateLayout() {},
      draw(params) {
        this.clear();
        this.addAxis();

        let domainFilter = (ann) =>
          ann.filter((a) => a.start < this.domain[1] && a.end > this.domain[0]);

        let domainWidth = this.domain[1] - this.domain[0];

        if (domainWidth < state.aliThresh) {
          let annotations = domainFilter(params.annotations);
          soda.sequence({
            chart: this,
            selector: "genome",
            annotations,
            row: 0,
          });
        }
      },
    });

    let segmentsRenderBlockLinks = (chart, block) => {
      if (block.link_data == undefined || block.link_data.length == 0)
        return undefined;
      let arc = soda.arc({
        chart: chart,
        annotations: block.link_data,
        strokeColor: "red",
        strokeWidth: 3,
        row: (d) => {
          return chart.layout.absRowToRow(d.a.row) - 1 + 14.5 / chart.rowHeight;
        },
        height: 12,
      });
      return arc;
    };

    let segments = new soda.Chart({
      ...chartConf,
      zoomable: true,
      rowHeight: 30,
      rowColors: ["whitesmoke", "white"],
      updateLayout(params) {
        // Must update the domain first inorder to properly layout the graph, by default this runs after this function but before draw
        // causing the graph to be out of sync...
        this.updateDomain(params);

        let query_to_row = new Map();
        //let fullRowCount = (state.showSegments)? Math.max(...params.historyBlocks.map((blk) => blk.query_id)) + 3: 0;

        let to_absolute_row = (blk) => (blk.row == 0 ? 0 : blk.query_id + 2);

        let visible_queries = params.historyBlocks
          .filter((blk) => {
            let start = blk.start;
            let end = blk.end;
            // Check if alignment is 'in bounds'...
            return end > this.domain[0] && start < this.domain[1];
          })
          .map(to_absolute_row)
          .sort((a, b) => a - b);

        let i = 0;
        for (const q_id of visible_queries) {
          if (!query_to_row.has(q_id)) {
            query_to_row.set(q_id, i);
            i += 1;
          }
        }

        let abs_row_to_row = (r) => {
          let r_floor = Math.floor(r);
          let new_r = query_to_row.get(r_floor) ?? -1;
          return new_r + (r % 1);
        };
        this.layout = {
          row: (d) => {
            let abs_row = to_absolute_row(d.a);
            return abs_row_to_row(abs_row);
          },
          toAbsRow: to_absolute_row,
          absRowToRow: abs_row_to_row,
          rowCount: state.showSegments ? query_to_row.size : 0,
        };
      },
      draw(params) {
        this.clear();
        if (this.showSegments) return;

        let y = (d) => this.rowHeight * this.layout.row(d) + 14;
        let x = (d) => this.xScale(d.a.start - 0.25);
        let width = (d) => this.xScale(d.a.end) - this.xScale(d.a.start - 0.5);
        let height = 13;

        let join_info = params.historyBlocks
          .filter((blk) => blk.segment != blk.join_to)
          .map((blk) => {
            return {
              row: this.layout.toAbsRow(blk),
              start: params.historySegments[blk.segment].end,
              end: params.historySegments[blk.join_to].start,
            };
          });

        function classColor(d) {
          if (d.a.row == 0) {
            return "#4d4d4d";
          }
          return classColorFromName(d.a.label);
        }

        // Color segments with alternating colors...
        soda.rectangle({
          chart: this,
          selector: "segments",
          annotations: params.historySegments,
          fillColor: (d) => (d.a.index % 2 ? "red" : "blue"),
          fillOpacity: 0.1,
          y: 0,
          x,
          width,
          height: this.viewportHeightPx,
        });

        // Render segment index and history count...
        soda.dynamicText({
          chart: this,
          selector: "segmentInfo",
          annotations: params.historySegments,
          row: 0,
          text: (d) => {
            return [`${d.a.index}: ${d.a.history_count}`, `${d.a.index}`];
          },
        });

        // Render the blocks...
        soda.rectangle({
          chart: this,
          selector: "blocks",
          annotations: params.historyBlocks,
          strokeColor: classColor,
          strokeWidth: 2,
          fillColor: "white",
          x,
          y,
          width,
          height,
        });

        // Display blocks that can join (how far ahead)...
        soda.arc({
          chart: this,
          selector: "blockJoins",
          annotations: join_info,
          row: (d) => {
            return this.layout.absRowToRow(d.a.row) - 1 + 14.5 / this.rowHeight;
          },
          height: 12,
        });

        // Display the name of the block...
        soda.dynamicText({
          chart: this,
          selector: "blocksLabel",
          annotations: params.historyBlocks,
          fontSize: 12,
          y: (d) => y(d) + 2,
          height,
          fontWeight: 500,
          fillColor: "black",
          text: (d) => {
            return [
              d.a.label,
              d.a.label.split("#")[0],
              d.a.label.charAt(0),
              "",
            ];
          },
        });

        for (let block of params.historyBlocks) {
          block.arc = segmentsRenderBlockLinks(this, block);
        }
      },
      postRender(params) {
        if (this.showSegments) return;

        let y = (d) => this.rowHeight * this.layout.row(d) + 14;
        let x = (d) => this.xScale(d.a.start - 0.25);
        let width = (d) => this.xScale(d.a.end) - this.xScale(d.a.start - 0.5);
        let height = 12;

        /*soda.hoverBehavior({
          chart: this,
          annotations: params.historyBlocks,
          // this function is evaluated when a glyph is moused over
          mouseover: (s, d) => s.style("stroke", "black"),
          // this function is evaluated when a glyph is no longer moused over
          mouseout: (s, d) => s.style("stroke", "none"),
          x,
          y,
          width,
          height,
        });*/

        soda.clickBehavior({
          chart: this,
          annotations: params.historyBlocks,
          click: (s, d) => {
            let processedLinks = d.a.link_data;

            if (processedLinks) {
              delete d.a.link_data;
              if (d.a.arc != undefined) {
                d.a.arc.remove();
                delete d.a.arc;
              }
            } else {
              let links = d.a.links;
              if (links == undefined || links.length == 0) return;

              let link_data = [];

              links.forEach((link) => {
                let other_segment = link.segment;
                if (other_segment == undefined) return;
                let [s1, s2] = [d.a.segment, other_segment].sort();

                link_data.push({
                  start: params.historySegments[s1].end,
                  end: params.historySegments[s2].start,
                  row: this.layout.toAbsRow(d.a),
                  weight: link.weight,
                });
              });

              d.a.link_data = link_data;
              d.a.arc = segmentsRenderBlockLinks(this, d.a);
            }
          },
        });

        soda.tooltip({
          chart: this,
          annotations: params.historySegments,
          row: 0,
          x,
          width,
          text: (d) => {
            return Object.entries(d.a)
              .map((val) => {
                let [k, v] = val;
                return `${k}: ${v}`;
              })
              .join("<br>\n");
          },
        });

        soda.tooltip({
          chart: this,
          annotations: params.historyBlocks,
          x,
          y,
          width,
          height,
          text: (d) => {
            return Object.entries(d.a)
              .map((val) => {
                let [k, v] = val;
                return `${k}: ${v}`;
              })
              .join("<br>\n");
          },
        });
      },
    });

    let alignments = new soda.Chart({
      ...chartConf,
      zoomable: true,
      rowHeight: 44,

      updateLayout(params) {
        let queryIds = [...new Set(params.proxy.map((a) => a.queryId))];
        let traceQueryIds = [];

        if (state.traceAtTop) {
          let allTrace = params.conclusiveTrace
            .concat(params.ambiguousTrace)
            .filter((a) => a.queryId != 0);

          allTrace.sort((a, b) => a.start - b.start);
          traceQueryIds = [...new Set(allTrace.map((a) => a.queryId))];
        }

        let remainingQueryIds = queryIds.filter(
          (id) => traceQueryIds.indexOf(id) === -1,
        );

        // the first two rows are reserved for the
        // skip state and the tandem repeat state
        let rowCount = 2;
        let dpRowToChartRow = new Map([[0, 0]]);
        params.tandemRepeats.forEach((r) => {
          dpRowToChartRow.set(r.row, 1);
        });

        let layoutFn = (id) => {
          let queryAssemblies = params.assemblies.filter(
            (a) => a.queryId == id,
          );
          let subLayout = soda.intervalGraphLayout(queryAssemblies);

          queryAssemblies.forEach((a) => {
            let subRow = subLayout.row({ a });
            dpRowToChartRow.set(a.row, rowCount + subRow);
          });

          rowCount += subLayout.rowCount;
        };

        traceQueryIds.forEach(layoutFn);
        remainingQueryIds.forEach(layoutFn);

        this.layout = {
          row: (d) => dpRowToChartRow.get(d.a.row) ?? -1,
          rowCount,
        };
      },

      draw(params) {
        let d3 = soda.internalD3;
        this.clear();
        // Doesn't work correctly....
        this.removeRowStripes();
        d3.select(this.highlightSelection.node().parentNode.parentNode).style(
          "background",
          `repeating-linear-gradient(to bottom, whitesmoke 0px, whitesmoke ${this.rowHeight}px, white  ${this.rowHeight}px, white ${this.rowHeight * 2}px)`,
        );

        let domainWidth = this.domain[1] - this.domain[0];

        let domainFilter = (ann) =>
          ann.filter((a) => a.start < this.domain[1] && a.end > this.domain[0]);

        let y = (d) => this.rowHeight * this.layout.row(d) + 14;
        let x = (d) => this.xScale(d.a.start - 0.25);
        let width = (d) => this.xScale(d.a.end) - this.xScale(d.a.start - 0.5);

        if (state.showInactive) {
          // inactive segments
          soda.rectangle({
            chart: this,
            selector: "inactive",
            annotations: params.inactiveSegments,
            fillColor: "red",
            fillOpacity: 0.1,
            y: 0,
            height: this.viewportHeightPx,
          });
        }

        // assemblies
        // NOTE: the y bug here seems to be when
        // the middle gap part of an assembly is
        // in the view, but the fragments aren't
        soda.rectangle({
          chart: this,
          selector: "assembly",
          annotations: params.assemblies.filter((a) => a.size > 1),
          y,
          width: (d) => this.xScale(d.a.end) - this.xScale(d.a.start + 1),
          height: 2,
          fillColor: (d) => {
            if (
              params.competedAssemblyRows[state.traceIteration].indexOf(
                d.a.row,
              ) >= 0
            ) {
              return options.competedColor;
            } else if (
              params.unresolvedAssemblyRows[state.traceIteration].indexOf(
                d.a.row,
              ) >= 0
            ) {
              return options.unresolvedColor;
            } else if (
              params.resolvedAssemblyRows[state.traceIteration].indexOf(
                d.a.row,
              ) >= 0
            ) {
              return options.resolvedColor;
            } else {
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
          fillColor: "white",
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

        // confidence segments
        soda.rectangle({
          chart: this,
          selector: "confidence-segments",
          annotations: params.confidenceSegments,
          x,
          y,
          width,
          height: 12,
          fillOpacity: (d) => d.a.conf,
          strokeColor: options.confidenceSegmentColor,
          strokeOpacity: 0.75,
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

          // tandem repeat labels
          soda.dynamicText({
            chart: this,
            selector: "tr-labels",
            annotations: params.tandemRepeats,
            fillColor: "black",
            text: (d) => [`tandem repeat(${d.a.period})`, "..."],
          });
        }

        if (domainWidth < state.aliThresh) {
          let annotations =
            params.alignmentConfidences != null
              ? domainFilter(params.alignmentConfidences)
              : [];
          annotations = annotations.map((val) => {
            let offset = val.start;
            let domain = this.domain;
            let start = Math.max(Math.floor(domain[0]), val.start);
            let end = Math.min(Math.ceil(domain[1]), val.end);
            return {
              id: val.id,
              row: val.row,
              values: val.values
                .slice(start - offset, end + 1 - offset)
                .map(Math.exp),
              scores: val.values.slice(start - offset, end + 1 - offset),
              norms: val.norms.slice(start - offset, end + 1 - offset),
              start: start - 0.5,
              end: end + 0.5,
            };
          });

          let heatmap_selection = soda.heatmap({
            chart: this,
            selector: "ali-seq-conf",
            annotations: annotations,
            y: (d) => y(d) + 12,
            height: 12,
          });

          heatmap_selection
            .selectAll(function () {
              return this.children;
            })
            .on("mouseover", function (data) {
              let element = this;
              let d3 = soda.internalD3;
              let bbox = element.getBoundingClientRect();
              let x = d3.event.clientX - bbox.left;
              let y = d3.event.clientY - bbox.top;

              let cell = Math.floor((x / bbox.width) * data.a.scores.length);

              d3.select(document.body)
                .select(".heatmap-hover-tooltip")
                .remove();

              d3.select(document.body)
                .select(".heatmap-highlight-tooltip")
                .remove();

              let cellX = Math.round(
                bbox.left + cell * (bbox.width / data.a.scores.length),
              );
              let cellY = bbox.top;
              let cellWidth = Math.round(bbox.width / data.a.scores.length);
              let cellHeight = Math.round(bbox.height);

              d3.select(document.body)
                .append("div")
                .attr("class", "heatmap-highlight-tooltip")
                .style("position", "fixed")
                .style("left", `${cellX}px`)
                .style("top", `${cellY}px`)
                .style("width", `${cellWidth}px`)
                .style("height", `${cellHeight}px`)
                .style("z-index", "1000")
                .style("pointer-events", "none")
                .style("background-color", "rgba(0, 255, 255, 0.4)");

              d3.select(document.body)
                .append("div")
                .attr("class", "heatmap-hover-tooltip")
                .style("position", "fixed")
                .style("z-index", "1000")
                .style("left", `${Math.round(cellX + cellWidth / 2)}px`)
                .style("top", `${cellY - 5}px`)
                .style("background-color", "lightblue")
                .style("border-radius", "3px")
                .style("transform", "translate(-50%, -100%)")
                .style("padding", "3px")
                .html(
                  `Log Score: ${data.a.scores[cell] + data.a.norms[cell]}<br>Normalized Log Score: ${data.a.scores[cell]}<br>Normalized Score: ${data.a.values[cell]}`,
                );
            })
            .on("mouseout", function () {
              let d3 = soda.internalD3;
              d3.select(document.body)
                .select(".heatmap-hover-tooltip")
                .remove();

              d3.select(document.body)
                .select(".heatmap-highlight-tooltip")
                .remove();
            });

          annotations = domainFilter(params.sequences);
          soda.sequence({
            chart: this,
            selector: "ali-seq",
            annotations,
            y: (d) => y(d) + 1,
            fillColor: (d) =>
              d.a.id[0] == "m" ? "green" : d.a.id[0] == "s" ? "orange" : "red",
          });

          annotations = domainFilter(params.gaps);
          soda.sequence({
            chart: this,
            selector: "ali-gaps",
            annotations,
            y: (d) => y(d) - 11,
            fillColor: "red",
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
          strokeColor: options.conclusiveColor,
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
          strokeColor: options.ambiguousColor,
          strokeWidth: 3,
        });

        soda.tooltip({
          annotations: params.ambiguousTrace.concat(params.conclusiveTrace),
          text: (d) => `confidence: ${d.a.conf}`,
        });

        soda.tooltip({
          annotations: params.confidenceSegments,
          text: (d) =>
            `${d.a.query}: ` +
            `${d.a.queryStart.toLocaleString()}..${d.a.queryEnd.toLocaleString()} / ${d.a.queryLength.toLocaleString()}` +
            `<br>chrom: ${d.a.start.toLocaleString()}..${d.a.end.toLocaleString()}` +
            `<br>strand: ${d.a.strand}` +
            `<br>confidence: ${d.a.conf}` +
            `<br>ali: ${d.a.ali_id}`,
        });
      },

      postZoom() {
        clearTimeout(timeoutId);
        timeoutId = window.setTimeout(() => {
          this.draw({
            ...this.renderParams,
            updateDomain: false,
          });

          genome.draw({
            ...genome.renderParams,
            updateDomain: false,
          });
        }, timeoutTime);
      },
    });

    alignments.render = function (params) {
      //this.resetTransform();

      let queryFilter = (a) => {
        // don't ever filter the skip state
        if (a.row == 0) {
          return true;
        }

        // this should prevent filtering TRs
        if (a.row > params.numQueries) {
          return true;
        }

        if (state.onlySelected) {
          if (a.end < params.start || a.start > params.end) {
            return false;
          }
        }

        if (state.onlyTrace) {
          if (state.traceRowsByIter[state.traceIteration].indexOf(a.row) < 0) {
            return false;
          }
        }

        if (state.regex != undefined) {
          let query = rowToQuery[a.row].toLowerCase();
          return state.regex.test(query);
        }

        if (rowToConf[a.row] <= state.confThresh) {
          return false;
        }

        return true;
      };

      let filteredParams = {
        ...params,
        assemblies: params.assemblies.filter(queryFilter),
        proxy: params.proxy.filter(queryFilter),
        sequences: params.sequences.filter(queryFilter),
        gaps: params.gaps.filter(queryFilter),
        ambiguousTrace:
          params.ambiguousTrace[state.traceIteration].filter(queryFilter),
        conclusiveTrace:
          params.conclusiveTrace[state.traceIteration].filter(queryFilter),
        inactiveSegments: params.inactiveSegments[state.traceIteration],
        confidenceSegments:
          params.confidenceSegments[state.traceIteration].filter(queryFilter),
        alignmentConfidences:
          params.alignmentConfidences != null
            ? params.alignmentConfidences.filter(queryFilter)
            : null,
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
    return {
      reference,
      referenceZoom,
      aurora,
      auroraZoom,
      genome,
      segments,
      alignments,
    };
  }

  function prepareAnn(ann) {
    let proxy = [];
    let aligned = [];
    let inner = [];
    let left = [];
    let right = [];
    let labels = [];
    let inversions = [];

    let inversionId = 0;

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

      let inversionLeft = -1;
      let inversionRight = -1;

      // aligned + inversions
      for (const [index, record] of group.aligned.entries()) {
        let tokens = record.split(",");
        let is_in_inversion =
          (index > 0 && group.strands[index] != group.strands[index - 1]) ||
          (index + 1 < group.strands.length &&
            group.strands[index] != group.strands[index + 1]);

        aligned.push({
          id: tokens[0],
          start: parseInt(tokens[1]),
          end: parseInt(tokens[2]),
          strand: group.strands[index],
          label: group.query,
          in_inversion: is_in_inversion,
        });

        if (is_in_inversion) {
          inversionLeft = inversionLeft < 0 ? index : inversionLeft;
          inversionRight = index;
        }
      }

      if (inversionLeft >= 0 && inversionRight >= 0) {
        inversions.push({
          id: "inv-" + inversionId,
          start: parseInt(group.aligned[inversionLeft].split(",")[1]),
          end: parseInt(group.aligned[inversionRight].split(",")[2]),
          label: group.query,
        });
        inversionId += 1;
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
      inversions,
    };
  }

  function prepareAli(ali) {
    let sequences = [];
    let proxy = [];
    let labelMap = new Map();
    let id = 0;
    let blank = "\u2000";
    let gaps = [];
    for (const a of ali) {
      let tokens = a.split(",");
      let green = tokens[0].replace(/ /g, blank);
      let orange = tokens[1].replace(/ /g, blank);
      let ali_gaps = [];
      if (tokens[2] != "") {
        ali_gaps = tokens[2].split("|").map((s) => {
          let [seq, start] = s.split(":");
          start = parseInt(start);
          return {
            start,
            end: start + seq.length,
            sequence: seq,
          };
        });
      }
      let start = parseInt(tokens[3]);
      let end = parseInt(tokens[4]);
      let query = tokens[5];
      let row = parseInt(tokens[6]);
      let queryId = parseInt(tokens[7]);
      let strand = tokens[8];

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

      let i = 0;
      for (let g of ali_gaps) {
        gaps.push({
          id: `g-${id}-${i++}`,
          row,
          queryId,
          ...g,
        });
      }

      id++;
    }
    return { sequences, gaps, proxy };
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
        let conf = parseFloat(tokens[4]);
        iterTrace.push({
          id: `trace-${iter}-${idx}`,
          traceIter: iter,
          start,
          end,
          queryId,
          row,
          conf,
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
      let period = parseInt(tokens[3]);
      let row = parseInt(tokens[4]);
      tandemRepeats.push({
        id: `tr-${row}`,
        start,
        end,
        consensus,
        period,
        row,
      });
    }
    return { tandemRepeats };
  }

  function prepareInactiveSegments(inactiveSegmentStrings) {
    let inactiveSegments = [];

    for (const [iter, iterStrings] of inactiveSegmentStrings.entries()) {
      let iterSegs = [];
      for (const [idx, seg] of iterStrings.entries()) {
        let tokens = seg.split(",");
        let start = parseInt(tokens[0]);
        let end = parseInt(tokens[1]);
        iterSegs.push({
          id: `ia-${iter}-${idx}`,
          start,
          end,
        });
      }
      inactiveSegments.push(iterSegs);
    }
    return { inactiveSegments };
  }

  function prepareConfidenceSegments(confidenceSegmentStrings) {
    let confidenceSegments = [];

    for (const [iter, iterStrings] of confidenceSegmentStrings.entries()) {
      let iterSegs = [];
      for (const [idx, seg] of iterStrings.entries()) {
        let tokens = seg.split(",");
        let start = parseInt(tokens[0]);
        let end = parseInt(tokens[1]);
        let row = parseInt(tokens[2]);
        let conf = parseFloat(tokens[3]);
        let queryStart = parseInt(tokens[4]);
        let queryEnd = parseInt(tokens[5]);
        let queryLength = parseInt(tokens[6]);
        let strand = tokens[7];
        let query = tokens[8];
        let ali_id = tokens[9];

        iterSegs.push({
          id: `cs-${iter}-${idx}`,
          start,
          end,
          row,
          conf,
          queryStart,
          queryEnd,
          queryLength,
          strand,
          query,
          ali_id,
        });
      }
      confidenceSegments.push(iterSegs);
    }
    return { confidenceSegments };
  }

  function prepareSegments(segmentStrings, targetStart) {
    let segments = [];

    for (const [index, seg] of segmentStrings.entries()) {
      let tokens = seg.split(",");
      let start = parseInt(tokens[0]) + targetStart;
      let end = parseInt(tokens[1]) + targetStart;
      let history_count = parseInt(tokens[2]);
      segments.push({
        id: `segment-${index}`,
        index,
        start,
        end,
        history_count,
      });
    }

    return segments;
  }

  function prepareBlocks(blockStrings, targetStart) {
    let blocks = [];

    for (const [index, blk] of blockStrings.entries()) {
      let tokens = blk.split(",");
      let segment = parseInt(tokens[0]);
      let block = parseInt(tokens[1]);
      let row = parseInt(tokens[2]);
      let query_id = parseInt(tokens[3]);
      let start = parseInt(tokens[4]) + targetStart;
      let end = parseInt(tokens[5]) + targetStart;
      let join_to = parseInt(tokens[6]);

      let confidence = parseFloat(tokens[7]);
      let alignment_score = parseFloat(tokens[8]);

      let label = tokens[9];
      let links = undefined;
      if (tokens[10].trim() != "") {
        links = tokens[10].split(";").map((v) => {
          let [segment, row, weight] = v.split(":");
          return {
            segment: parseInt(segment),
            row: parseInt(row),
            weight: parseFloat(weight),
          };
        });
      } else {
        links = [];
      }

      blocks.push({
        id: `block-${index}`,
        start,
        end,
        segment,
        block,
        row,
        query_id,
        join_to,
        confidence,
        alignment_score,
        label,
        links,
      });
    }

    return blocks;
  }

  function isLittleEndian() {
    let arr = new Uint32Array([0x11223344]);
    let view = new Uint8Array(arr.buffer);
    return view[0] == 0x44;
  }

  function base64ToFloats(data) {
    let dataString = atob(data);
    let intView = new Uint8Array(dataString.length);
    for (let i = 0; i < dataString.length; i++)
      intView[i] = dataString.charCodeAt(i);
    if (!isLittleEndian()) {
      for (let i = 0; i < intView.length; i += 4) {
        for (let j = 0; j < 4; j++) {
          let tmp = intView[i + (3 - j)];
          intView[i + (3 - j)] = intView[i + j];
          intView[i + j] = tmp;
        }
      }
    }
    return new Float32Array(intView.buffer);
  }

  function prepareAlignmentConfidences(alScores) {
    if (alScores == null) return alScores;

    let alignment_scores = [];
    let i = 0;

    for (const entry of alScores) {
      let tokens = entry.split(",");

      let start = parseInt(tokens[0]);
      let end = parseInt(tokens[1]);
      let values = base64ToFloats(tokens[2]);
      let norms = base64ToFloats(tokens[3]);

      alignment_scores.push({
        start,
        end,
        values,
        norms,
        id: `${i}`,
        row: i,
        average: values.reduce((acc, v) => acc + v, 0) / values.length,
      });

      i++;
    }

    return alignment_scores;
  }

  function prepareData() {
    let coords = {
      start: data.targetStart - LABEL_WIDTH,
      end: data.targetEnd + LABEL_WIDTH,
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

    let genomeAnn = [];
    let cnt = 0;
    for (let i = 0; i < data.targetSeq.length; i += 1000) {
      let start = data.targetStart + i;
      let end = Math.min(data.targetEnd + 1, start + 1000);
      let split = data.targetSeq.slice(i, i + 1000);
      genomeAnn.push({
        id: `gseq-${cnt++}`,
        start,
        end,
        sequence: split,
      });
    }

    let genome = {
      ...coords,
      annotations: genomeAnn,
    };

    let alignments = {
      ...coords,
      numQueries: data.numQueries,
      ...prepareAli(data.alignmentStrings),
      ...prepareAssemblies(data.assemblyStrings),
      ...prepareTandemRepeats(data.tandemRepeatStrings),
      ambiguousTrace: prepareTrace(
        data.ambiguousTraceStrings,
        data.targetStart,
      ),
      conclusiveTrace: prepareTrace(
        data.conclusiveTraceStrings,
        data.targetStart,
      ),
      resolvedAssemblyRows: data.resolvedAssemblyRows,
      unresolvedAssemblyRows: data.unresolvedAssemblyRows,
      competedAssemblyRows: data.competedAssemblyRows,
      ...prepareInactiveSegments(data.inactiveSegmentStrings),
      ...prepareConfidenceSegments(data.confidenceSegmentStrings),
      historySegments: prepareSegments(data.historySegments, data.targetStart),
      historyBlocks: prepareBlocks(data.historyBlocks, data.targetStart),
      alignmentConfidences: prepareAlignmentConfidences(
        data.alignmentConfidences,
      ),
    };

    alignments.proxy.forEach((a) => {
      rowToQuery[a.row] = a.query;
      rowToConf[a.row] = 0.0;
    });

    alignments.confidenceSegments[0].forEach((a) => {
      rowToConf[a.row] = Math.max(rowToConf[a.row], a.conf);
    });

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

    charts.referenceZoom.render({
      ...params.reference,
      ...coords,
    });

    charts.auroraZoom.render({
      ...params.aurora,
      ...coords,
    });

    charts.segments.render({
      ...params.alignments,
      labels: params.aurora.labels,
      ...coords,
    });

    charts.alignments.render({
      updateDomain,
      ...params.alignments,
      ...coords,
    });

    // Remove generated svg elements...
    for (let emptySvg of document.querySelectorAll("body > svg")) {
      emptySvg.remove();
    }
  }

  function render() {
    charts.reference.render(params.reference);
    charts.aurora.render(params.aurora);
    initializeBrush();

    let zoomSync = new soda.ZoomSyncer();
    zoomSync.add([
      charts.auroraZoom,
      charts.referenceZoom,
      charts.genome,
      charts.segments,
      charts.alignments,
    ]);
  }
}

document.addEventListener("DOMContentLoaded", function() {
  bootstrap();
});
