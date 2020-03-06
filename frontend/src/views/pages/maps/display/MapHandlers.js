import React from 'react';

class MapHandlers extends React.Component {
  constructor(props) {
    super(props);

    this.state = {
      hoverMol: null
    }
  }

  handleMolsSelect = (mols, points) => {
    if (this.props.onMolsSelect) {
      this.props.onMolsSelect(mols);
    }
    if (this.props.onPointsSelect) {
      this.props.onPointsSelect(points);
    }
  };

  handleMolHover = (mol, point) => {
    this.setState({
      hoverMol : mol
    })
  };

  render() {
    const Comp = this.props.component;
    return (<Comp
      {...this.props}
      {...this.state}
      onMolsSelect={this.handleMolsSelect}
      onMolHover={this.handleMolHover}
    />)
  }
}

export default MapHandlers;