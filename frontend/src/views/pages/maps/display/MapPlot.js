import React from 'react';
import Plot from 'react-plotly.js';
// import { Popover, PopoverBody, PopoverHeader } from 'reactstrap';
// import { MoleculeImage } from '../../../../genui';

// function MoleculePopover(props) {
//   return props.mol ? (
//     <Popover placement="top" isOpen={props.open} target="mol-popover">
//       <PopoverHeader>Closest Molecule</PopoverHeader>
//       <PopoverBody>
//         <MoleculeImage mol={props.mol}/>
//       </PopoverBody>
//     </Popover>
//   ) : null
// }

class MapPlot extends React.Component {

  constructor(props) {
    super(props);

    const map = this.props.map;

    this.layout = {
      title: map.name,
      xaxis: this.initAxis(map, 'x'),
      yaxis: this.initAxis(map, 'y'),
      showlegend: true,
      legend: {
        x: 1,
        xanchor: 'right',
        y: 1
      },
      font: {size: 18},
      autosize: true,
      dragmode: 'lasso',
      hovermode: 'closest',
      datarevision: 0,
    };

    this.config = {
      responsive: false,
      displaylogo: false,
      displayModeBar: true
    };

    this.state = {
      traces : this.initTraces(map),
      layout: this.layout,
      revision: 0,
      popover: {
        open: false
        , x: 1
        , y: 1
        , mol : null
      }
    }
  }

  initAxis = (map, axis) => {
    const algorithm = map.trainingStrategy.algorithm;
    return {
      title: {
        text: `${algorithm.name}-${axis}`,
        // font: {
        //   family: 'Courier New, monospace',
        //   size: 18,
        //   color: '#7f7f7f'
        // }
      }
    }
  };

  initTraces = (map) => {
    return map.molsets.map(molset => ({
      x: [],
      y: [],
      mode: 'markers',
      type: 'scatter',
      name: molset.name,
      customdata: [],
      marker: {
        color: this.props.molsetsToColor[molset.id],
      }
    }))
  };

  shouldComponentUpdate(nextProps, nextState, nextContext) {
    if (this.plotlyRef && this.plotlyRef.resizeHandler) {
      this.plotlyRef.resizeHandler();
    }

    return true;
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    const points = this.props.points;
    const map = this.props.map;
    const molsets = map.molsets;

    if (prevProps.map.id !== map.id) {
      prevState.layout = Object.assign(prevState.layout, {
        title: map.name,
        xaxis: this.initAxis(map, 'x'),
        yaxis: this.initAxis(map, 'y'),
      });
      prevState.traces = this.initTraces(map);
      this.setState(prevState);
      return
    }

    Object.entries(points).forEach(entry => {
      const val = entry[1];
      const key = Number(entry[0]);
      const traceIdx = molsets.findIndex(molset => molset.id === key);
      const prevTrace = prevState.traces[traceIdx];
      const prevLength = prevTrace.x.length;
      const newTrace = this.updateTrace(prevTrace, val);
      if (newTrace.x.length !== prevLength) {
        prevState.traces[traceIdx] = newTrace;
        this.setState(prev => {
          prevState.layout.datarevision = prev.layout.datarevision + 1;
          prevState.revision = prev.revision + 1;
          return prevState
        })
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
    if (eventData && this.props.onMolsSelect) {
      this.props.onMolsSelect(
        eventData.points.map(point => point.customdata.molecule),
        eventData.points,
      );

      if (this.props.setSelectedMolsInMapRevision) {
        this.props.setSelectedMolsInMapRevision(this.props.selectedMolsRevision + 1);
      }

    }
  };

  handleDeselect = () => {
    if (this.props.onDeselect) {
      this.props.onDeselect();
      if (this.props.setSelectedMolsInMapRevision) {
        this.props.setSelectedMolsInMapRevision(this.props.selectedMolsRevision + 1);
      }
    }
  };

  handleHover = (eventData) => {
    const point = eventData.points[0];
    const data = {
      open : true,
      x : eventData.event.clientX,
      y : eventData.event.clientY,
      mol: point.customdata.molecule,
    };
    this.setState({popover : data}
    );
    this.props.onMolHover(data.mol, point)
  };

  render() {
    const traces = this.state.traces;
    return (
      <div className="genui-map-plot">
        <Plot
          ref={plotlyRef => {
            this.plotlyRef = plotlyRef;
          }}
          data={traces}
          revision={this.state.revision}
          layout={this.layout}
          config={this.config}
          useResizeHandler={true}
          style={{width: "100%", height: "100%"}}
          onSelected={this.handleSelect}
          onHover={this.handleHover}
          onDeselect={this.handleDeselect}
          onSelecting={() => {
            this.setState({popover : {
                open : false
              }});
          }}
        />
        {/*<div*/}
        {/*  id="mol-popover"*/}
        {/*>*/}
        {/*  <MoleculePopover {...this.state.popover} />*/}
        {/*</div>*/}
      </div>
    )
  }
}

export default MapPlot;