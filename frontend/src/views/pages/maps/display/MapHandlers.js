import React from 'react';

class MapHandlers extends React.Component {
  constructor(props) {
    super(props);

    this.state = {
      selectedMols : [],
      selectedPoints : [],
      hoverMol: null
    }
  }

  handleMolsSelect = (mols, points) => {
    this.setState({
      selectedMols : mols,
      selectedPoints: points,
    })
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