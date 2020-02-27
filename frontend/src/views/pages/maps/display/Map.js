import React from 'react';
import Plot from 'react-plotly.js';
import './plot-styles.css';

export default function Map(props) {
  let trace1 = {
    x: [1, 2, 3, 4],
    y: [10, 15, 13, 17],
    mode: 'markers',
    type: 'scatter'
  };

  let trace2 = {
    x: [2, 3, 4, 5],
    y: [16, 5, 11, 9],
    mode: 'lines',
    type: 'scatter'
  };

  let trace3 = {
    x: [1, 2, 3, 4],
    y: [12, 9, 15, 12],
    mode: 'lines+markers',
    type: 'scatter'
  };

  let data = [trace1, trace2, trace3];

  let layout = {
    title: 'Responsive to window\'s size!',
    font: {size: 18},
    autosize: true
  };

  let config = {responsive: true};

  return (
    <div className="genui-map-plot">
      <Plot
        data={data}
        layout={layout}
        config={config}
        useResizeHandler={true}
        style={{width: "100%", height: "100%"}}
      />
    </div>
  )
}