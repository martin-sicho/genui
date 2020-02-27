import React from 'react';
import Plot from 'react-plotly.js';
import './plot-styles.css';
import withUnmounted from '@ishawnwang/withunmounted';

function MapPlot(props) {

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
    title: 'Responsive',
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

class Map extends React.Component {
  abort = new AbortController();
  hasUnmounted = false;

  constructor(props) {
    super(props);

    this.map = this.props.map;
    this.molsets = this.map.molsets;
    this.pointsUrl = new URL(`${this.map.id}/points/`, this.props.apiUrls.mapsRoot);

    this.state = {
      points : {},
      lastPage : null,
      nextPage : this.pointsUrl,
    }
  }

  componentDidMount() {
    this.fetchData(this.state.nextPage);
  }

  componentWillUnmount() {
    this.abort.abort();
  }

  fetchData = (page) => {
    fetch(page, {signal : this.abort.signal})
      .then(response => this.props.handleResponseErrors(response, `Cannot fetch page: ${page.toString()}`))
      .then(data => {
        const pagePoints = {};
        data.results.forEach(
          result => {
            const providers = [];
            this.molsets.forEach(molset => {
              if (result.molecule.providers.includes(molset.id)) {
                providers.push(molset);
              }
            });
            providers.forEach(provider => {
              if (!pagePoints.hasOwnProperty(provider.id)) {
                pagePoints[provider.id] = [];
              }
              pagePoints[provider.id].push(result)
            })
          }
        );

        if (this.hasUnmounted) {
          return
        }

        let nextPage = null;
        if (data.next) {
          nextPage = new URL(data.next);
          this.fetchData(new URL(data.next));
        }

        this.setState(prevState => {
          const points = prevState.points;
          Object.keys(pagePoints).forEach(providerID => {
            if (points[providerID]) {
              pagePoints[providerID].forEach(result =>
                points[providerID].push(result)
              );
            } else {
              points[providerID] = pagePoints[providerID];
            }
          });

          return {
            points : points,
            lastPage: page,
            nextPage: nextPage,
          }
        })
      })
      .catch(e => console.log(e))
  };

  render() {
    return (
      <MapPlot
        {...this.props}
        molsets={this.molsets}
        points={this.state.points}
      />
    )
  }
}

export default withUnmounted(Map);