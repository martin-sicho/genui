import React from 'react';
import { Col, Row } from 'reactstrap';
import { ActivitySetList, MoleculeActivityDetail, MoleculeMetadata, MoleculeImage } from '../../..';
import SimplePaginator from '../../SimplePaginator';

export function CompoundListItem(props) {
  const mol = props.mol;
  const showData = typeof(props.showInfo) === 'boolean' ? props.showInfo : true;
  const showActivities = typeof(props.showActivities) === 'boolean' ? props.showActivities : true;
  const sm_cols = [3, 3, 6];
  const md_cols = [3, 3, 6];
  if (!showData) {
    sm_cols[0] = sm_cols[0] + 1;
    md_cols[0] = md_cols[0] + 1;
    sm_cols[2] = sm_cols[2] + 2;
    md_cols[2] = md_cols[2] + 2;
  }
  if (!showActivities) {
    sm_cols[0] = sm_cols[0] + sm_cols[2] / 2;
    md_cols[0] = md_cols[0] + md_cols[2] / 2;
    sm_cols[1] = sm_cols[1] + sm_cols[2] / 2;
    md_cols[1] = md_cols[1] + md_cols[2] / 2;
  }

  return (
    <Row>
      <Col md={md_cols[0]} sm={sm_cols[0]}>
        <MoleculeImage mol={mol}/>
      </Col>
      {
        showData ? (
          <Col md={md_cols[1]} sm={sm_cols[1]}>
            <MoleculeMetadata {...props} mol={mol}/>
          </Col>
        ) : null
      }
      {
        showActivities ? (
          <Col md={md_cols[2]} sm={sm_cols[2]}>
            <MoleculeActivityDetail
              {...props}
              mol={mol}
              component={ActivitySetList}
            />
          </Col>
        ) : null
      }
    </Row>
  )
}

export function CompoundListPageItem(props) {
  return <CompoundListItem {...props} mol={props.pageItem}/>
}

export default function CompoundList(props) {
  const mols = props.mols;

  if (props.paginate) {
    return (
      <SimplePaginator {...props} items={mols} component={CompoundListPageItem}/>
    )
  } else {
      return (
        <React.Fragment>
          {
            mols.map(mol => (
              <CompoundListItem {...props} key={mol.id} mol={mol}/>
            ))
          }
        </React.Fragment>
      )
  }
}