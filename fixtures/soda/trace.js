import * as soda from "https://cdn.skypack.dev/@sodaviz/soda@0.13.1";

export function run(data) {
  let chart = new soda.Chart({
    selector: "div.container",
    zoomable: true,
    resizable: true,
    rowHeight: 20,
    rowColors: ["gainsboro", "white"],
    draw(params) {
      this.addAxis();

      soda.rectangle({
        chart: this,
        annotations: params.alignments,
        fillColor: "red",
        fillOpacity: 0.2,
        row: (d) => d.a.row,
      });

      soda.tooltip({
        chart: this,
        annotations: params.alignments,
        text: (d) => d.a.name,
      });

      soda.line({
        chart: this,
        annotations: params.lines,
        row: (d) => d.a.row,
      });

      soda.line({
        chart: this,
        annotations: params.spawns,
        x1: (d) => d.c.xScale(d.a.col1),
        y1: (d) => d.c.rowHeight * d.a.row1 + 10,
        x2: (d) => d.c.xScale(d.a.col2),
        y2: (d) => d.c.rowHeight * d.a.row2 + 10,
      });
    },
  });

  let lines = data.lines.map((s) => {
    let tokens = s.split(",");
    return {
      id: `l-${tokens[0]}`,
      start: parseInt(tokens[1]),
      end: parseInt(tokens[2]),
      row: parseInt(tokens[3]),
    };
  });

  let spawns = data.spawns.map((s) => {
    let tokens = s.split(",");
    return {
      id: `s-${tokens[0]}`,
      row1: parseInt(tokens[1]),
      col1: parseInt(tokens[2]),
      row2: parseInt(tokens[3]),
      col2: parseInt(tokens[4]),
    };
  });

  let alignments = data.alignments.map((a) => {
    let tokens = a.split(",");
    return {
      id: `a-${tokens[0]}`,
      name: tokens[1],
      start: parseInt(tokens[2]),
      end: parseInt(tokens[3]),
      row: parseInt(tokens[0]) + 1,
    };
  });

  chart.render({
    lines,
    spawns,
    alignments,
    start: data.start,
    end: data.end,
    rowCount: data.rowCount,
  });
}
