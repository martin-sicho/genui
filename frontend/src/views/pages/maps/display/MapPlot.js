import React from 'react';
import Plot from 'react-plotly.js';

class MapPlot extends React.Component {

  constructor(props) {
    super(props);

    this.map = this.props.map;
    this.molsets = this.map.molsets;

    const algorithm = this.map.trainingStrategy.algorithm;
    this.layout = {
      title: this.map.name,
      font: {size: 18},
      autosize: true,
      dragmode: 'lasso',
      xaxis: {
        title: {
          text: `${algorithm.name}-x`,
          // font: {
          //   family: 'Courier New, monospace',
          //   size: 18,
          //   color: '#7f7f7f'
          // }
        },
      },
      yaxis: {
        title: {
          text: `${algorithm.name}-y`,
          // font: {
          //   family: 'Courier New, monospace',
          //   size: 18,
          //   color: '#7f7f7f'
          // }
        }
      }
    };

    this.config = {
      responsive: true,
      displaylogo: false,
      displayModeBar: true
    };

    this.state = {
      traces : this.molsets.map(molset => ({
        x: [],
        y: [],
        mode: 'markers',
        type: 'scatter',
        name: molset.name,
        customdata: [],
      })),
      layout: this.layout,
      revision: 0
    }
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    const points = this.props.points;
    const molsets = this.molsets;

    Object.entries(points).forEach(entry => {
      const val = entry[1];
      const key = Number(entry[0]);
      const traceIdx = molsets.findIndex(molset => molset.id === key);
      const prevTrace = prevState.traces[traceIdx];
      const prevLength = prevTrace.x.length;
      const newTrace = this.updateTrace(prevTrace, val);
      if (newTrace.x.length !== prevLength) {
        prevState.traces[traceIdx] = newTrace;
        prevState.revision++;
        prevState.layout.datarevision++;
        this.setState(prevState)
      }
    });
  }

  updateTrace = (prevTrace, newPoints) => {
    const nextIdx = prevTrace.x.length;
    for (let idx = nextIdx; idx < newPoints.length;idx++) {
      const newPoint = newPoints[idx];
      prevTrace.x.push(newPoint.x);
      prevTrace.y.push(newPoint.y);
      prevTrace.customdata.push(newPoint)
    }
    return prevTrace;
  };

  handleSelect = (eventData) => {
    if (this.props.onMolsSelect) {
      this.props.onMolsSelect(
        eventData.points.map(point => point.customdata.molecule),
        eventData.points,
      );
    }
  };

  render() {
    const traces = this.state.traces;
    return (
      <div className="genui-map-plot">
        <Plot
          data={traces}
          revision={this.state.revision}
          layout={this.layout}
          config={this.config}
          useResizeHandler={true}
          style={{width: "100%", height: "100%"}}
          onSelected={this.handleSelect}
        />
      </div>
    )
  }
}

export default MapPlot;