import * as soda from "https://cdn.skypack.dev/@sodaviz/soda@0.13.1";
import * as d3 from "https://esm.run/d3-scale@4";

let colors = [
  "#e6194b",
  "#3cb44b",
  "#ffe119",
  "#4363d8",
  "#f58231",
  "#911eb4",
  "#f032e6",
  "#bcf60c",
  "#fabebe",
  "#008080",
  "#e6beff",
  "#9a6324",
  "#fffac8",
  "#800000",
  "#aaffc3",
  "#808000",
  "#ffd8b1",
  "#000075",
  "#808080",
  "#000000",
];

let colorScale = d3.scaleLinear().domain([0.0, 1.0]).range(["white", "black"]);

export function run(data) {
  let targetChart = new soda.Chart({
    selector: "div.container",
    zoomable: true,
    resizable: true,
    rowHeight: 30,
    draw(params) {
      this.addAxis();

      soda.rectangle({
        chart: this,
        selector: "target",
        annotations: params.annotations,
        strokeWidth: 2,
        strokeColor: (d) => colors[parseInt(d.a.id) % colors.length],
        fillColor: (d) => colorScale(d.a.conf),
        row: 1,
      });

      soda.arc({
        chart: this,
        selector: "links",
        annotations: params.links,
        row: 0,
      });

      soda.tooltip({
        chart: this,
        annotations: params.annotations,
        text: (d) => `${d.a.conf.toFixed(4)}`,
      });
    },
  });

  let consensusChart = new soda.Chart({
    selector: "div.container",
    zoomable: true,
    resizable: true,
    rowHeight: 10,
    draw(params) {
      this.addAxis();

      soda.rectangle({
        chart: this,
        strokeWidth: 2,
        strokeColor: (d) => colors[parseInt(d.a.id) % colors.length],
        fillColor: (d) => colors[parseInt(d.a.id) % colors.length],
        selector: "consensus",
        annotations: params.annotations,
      });

      soda.hoverBehavior({
        annotations: params.annotations,
        mouseover: (s, d) => {
          let glyphs = soda.queryGlyphMap({ annotations: [d.a] });
          glyphs.forEach((g) => g.style("stroke", "cyan"));
        },
        mouseout: (s, d) => {
          let glyphs = soda.queryGlyphMap({ annotations: [d.a] });
          glyphs.forEach((g) =>
            g.style("stroke", colors[parseInt(d.a.id) % colors.length]),
          );
        },
      });
    },
  });

  let assemblyTargetConf = {
    selector: "div.container",
    zoomable: true,
    resizable: true,
    rowHeight: 30,
    divOutline: "1px solid black",
    draw(params) {
      soda.rectangle({
        chart: this,
        strokeWidth: 2,
        strokeColor: (d) => colors[parseInt(d.a.id) % colors.length],
        fillColor: (d) => colorScale(d.a.conf),
        selector: "target",
        annotations: params.annotations,
      });
    },
  };

  let assemblyConsensusConf = {
    selector: "div.container",
    zoomable: true,
    resizable: true,
    rowHeight: 10,
    divOutline: "1px solid black",
    draw(params) {
      soda.rectangle({
        chart: this,
        strokeWidth: 2,
        strokeColor: (d) => colors[parseInt(d.a.id) % colors.length],
        fillColor: (d) => colors[parseInt(d.a.id) % colors.length],
        selector: "consensus",
        annotations: params.annotations,
      });
    },
  };

  let targetParams = {
    annotations: data.targetAliStrings.map((r) => {
      let tokens = r.split(",")
      return {
        id: tokens[0],
        start: parseFloat(tokens[1]),
        end: parseFloat(tokens[2]),
        conf: parseFloat(tokens[3]),
      }
    }),
    links: data.links.map((r) => {
      let tokens = r.split(",")
      return {
        id: tokens[0],
        start: parseFloat(tokens[1]),
        end: parseFloat(tokens[2]),
      }
    }),
    rowCount: 2,
  };

  let consensusParams = {
    annotations: data.consensusAliStrings.map((r) => {
      let tokens = r.split(",")
      return {
        id: tokens[0],
        start: parseFloat(tokens[1]),
        end: parseFloat(tokens[2]),
        conf: parseFloat(tokens[3]),
      }
    }),
  };

  let targetAssemblyCharts = [];
  let consensusAssemblyCharts = [];
  let numAssemblies = data.targetAssemblyStrings.length;
  
  console.log(data);
  for (let i = 0; i < numAssemblies; i++) {
    let targetChart = new soda.Chart(assemblyTargetConf);
    targetChart.render({
      start: data.targetStart,
      end: data.targetEnd,
      annotations: data.targetAssemblyStrings[i].map((r) => {
        let tokens = r.split(",")
        return {
          id: tokens[0],
          start: parseFloat(tokens[1]),
          end: parseFloat(tokens[2]),
          conf: parseFloat(tokens[3]),
        }
      }),
    });
    targetAssemblyCharts.push(targetChart);

    let consensusChart = new soda.Chart(assemblyConsensusConf);
    consensusChart.render({
      start: data.consensusStart,
      end: data.consensusEnd,
      annotations: data.consensusAssemblyStrings[i].map((r) => {
        let tokens = r.split(",")
        return {
          id: tokens[0],
          start: parseFloat(tokens[1]),
          end: parseFloat(tokens[2]),
          conf: parseFloat(tokens[3]),
        }
      }),
    });
    consensusAssemblyCharts.push(consensusChart);
  }
  
  targetChart.render(targetParams);
  consensusChart.render(consensusParams);

  let zsTarget = new soda.ZoomSyncer();
  zsTarget.addChart(targetChart);
  zsTarget.add(targetAssemblyCharts);

  let zsConsensus = new soda.ZoomSyncer();
  zsConsensus.addChart(consensusChart);
  zsConsensus.add(consensusAssemblyCharts);
}
