import React from 'react';
import './plot-styles.css';
import withUnmounted from '@ishawnwang/withunmounted';
import MapPlot from './MapPlot';

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
        points={this.state.points}
      />
    )
  }
}

export default withUnmounted(Map);